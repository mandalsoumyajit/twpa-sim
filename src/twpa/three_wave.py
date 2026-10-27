"""Flux-biased rf-SQUID/SNAIL chains: finite-circuit 3WM and Gaussian noise.

Node phases are departures from the selected zero-current equilibrium.
The exact potential is projected in frequency space, including DC rectification.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve
from .models import positive
from .jpa_circuit import PHIBAR, PHI0, Tone, Transient, voltage_sources
from .jtwpa import (Chain, ResonantLoad, _basis, _linear_system, _newton,
                     _jacobian, _real_vector, _complex_vector)
from .jtwpa import bloch as _bloch, linear_response as _linear_response


@dataclass(frozen=True)
class RFSQUID:
    critical_current: float
    loop_inductance: float
    junction_capacitance: float
    junction_resistance: float = 1e9
    flux_bias: float = .3  # Phi/Phi0.

    def __post_init__(self):
        for key in ('critical_current', 'loop_inductance', 'junction_capacitance', 'junction_resistance'):
            positive(getattr(self, key), key)
        if not np.isfinite(self.flux_bias):
            raise ValueError('flux_bias must be finite')
        if self.beta >= 1:
            raise ValueError('rf-SQUID requires beta_L < 1 for its unique supported equilibrium')

    @property
    def beta(self):
        return self.loop_inductance*self.critical_current/PHIBAR

    @property
    def equilibrium(self):
        ext = 2*np.pi*self.flux_bias
        return brentq(lambda x: x+self.beta*np.sin(x)-ext, ext-self.beta-1e-12, ext+self.beta+1e-12)

    def current(self, phase, derivative=0):
        if derivative == 0:
            return self.critical_current*(np.sin(phase)+(phase-2*np.pi*self.flux_bias)/self.beta)
        if derivative == 1:
            return self.critical_current*(np.cos(phase)+1/self.beta)
        return self.critical_current*np.sin(phase+derivative*np.pi/2)

    @property
    def capacitance(self):
        return self.junction_capacitance

    @property
    def resistance(self):
        return self.junction_resistance


@dataclass(frozen=True)
class SNAIL:
    critical_current: float  # Each large junction.
    asymmetry: float
    junctions: int
    junction_capacitance: float  # Each large junction; small scales with asymmetry.
    junction_resistance: float = 1e9  # Small R = large R / asymmetry.
    flux_bias: float = .35

    def __post_init__(self):
        for key in ('critical_current', 'junction_capacitance', 'junction_resistance'):
            positive(getattr(self, key), key)
        if not isinstance(self.junctions, int) or self.junctions < 2:
            raise ValueError('SNAIL needs at least two large junctions')
        if not np.isfinite(self.asymmetry) or not 0 < self.asymmetry < 1/self.junctions:
            raise ValueError('initial SNAIL model requires 0 < asymmetry < 1/junctions')
        if not np.isfinite(self.flux_bias):
            raise ValueError('flux_bias must be finite')
        self.equilibrium  # Reject ambiguous/unstable branches immediately.

    @property
    def equilibrium(self):
        ext, n = 2*np.pi*self.flux_bias, self.junctions
        grid = ext+np.linspace(-n*np.pi/2, n*np.pi/2, 513)
        values = self.current(grid)
        roots = []
        for lo, hi, a, b in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if a == 0 or a*b < 0:
                x = lo if a == 0 else brentq(self.current, lo, hi, xtol=1e-14)
                if self.current(x, 1)>0 and not any(abs(x-y)<1e-9 for y in roots):
                    roots.append(x)
        if len(roots) != 1:
            raise ValueError('SNAIL equilibrium is ambiguous; select a different supported flux/asymmetry')
        return roots[0]

    def current(self, phase, derivative=0):
        n, ext = self.junctions, 2*np.pi*self.flux_bias
        return self.critical_current*(self.asymmetry*np.sin(phase+derivative*np.pi/2)
              +np.sin((phase-ext)/n+derivative*np.pi/2)/n**derivative)

    @property
    def capacitance(self):
        return (self.asymmetry+1/self.junctions)*self.junction_capacitance

    @property
    def resistance(self):
        return self.junction_resistance/(self.asymmetry+1/self.junctions)


def coefficients(element):
    """Potential derivatives c2,c3,c4 in units of large/single-junction EJ.

    U/EJ = constant + c2*x^2/2 + c3*x^3/6 + c4*x^4/24 + ... .
    """
    eq = element.equilibrium
    return np.array([element.current(eq, order)/element.critical_current for order in (1, 2, 3)])


@dataclass(frozen=True)
class ThreeWaveChain:
    element: RFSQUID | SNAIL
    ground_capacitance: float
    cells: int
    port_impedance: float = 50.
    ground_resistance: float | None = None
    cell_length: float = 10e-6
    resonator: ResonantLoad | None = None

    def __post_init__(self):
        if not isinstance(self.element, (RFSQUID, SNAIL)):
            raise TypeError('element must be RFSQUID or SNAIL')
        self.linear_chain

    @property
    def linear_chain(self):
        e = self.element
        return Chain(e.current(e.equilibrium, 1), self.ground_capacitance, e.capacitance,
                     self.cells, self.port_impedance, e.resistance, self.ground_resistance,
                     self.cell_length, self.resonator)


def bloch(chain, frequency):
    return _bloch(chain.linear_chain, frequency)


def linear_response(chain, frequency):
    return _linear_response(chain.linear_chain, frequency)


class _Potential:
    def __init__(self, element):
        self.element = element
        self.equilibrium = element.equilibrium
        self.scale = element.current(self.equilibrium, 1)

    def dc(self, waveform):
        offset = np.zeros(waveform.shape[1])
        for _ in range(30):
            phase = self.equilibrium+waveform+offset
            current = self.element.current(phase)/self.scale
            curvature = self.element.current(phase, 1)/self.scale
            mean_curvature = np.mean(curvature, axis=0)
            if np.any(mean_curvature <= 1e-5):
                raise ValueError('pump leaves the supported positive-curvature DC branch')
            error = np.mean(current, axis=0)
            if max(abs(error)) < 2e-13:
                return offset, current, curvature
            offset -= np.clip(error/mean_curvature, -.3, .3)
        raise RuntimeError('pump-induced DC phase did not converge')

    def project(self, branch, indices, quadrature, derivative):
        basis = _basis(tuple(indices), quadrature)
        waveform = np.real(basis.T@branch)
        _, current, curvature = self.dc(waveform)
        projection = 2*basis.conj()/basis.shape[1]
        if not derivative:
            return branch-projection@current
        result = []
        for column in range(2*len(branch)):
            local = basis[column % len(branch)].real if column<len(branch) else -basis[column % len(branch)].imag
            # Static rectification is eliminated at zero mean current. This
            # Schur complement retains a symmetric conservative Hessian.
            dc_derivative = -np.mean(curvature*local[:, None], axis=0)/np.mean(curvature, axis=0)
            value = -projection@(curvature*(local[:, None]+dc_derivative))
            value[column % len(branch)] += 1 if column<len(branch) else 1j
            result.append(value)
        return result


def _state(chain, pump_frequency, pump_power, signal_frequency, *, harmonics=3,
           sidebands=1, quadrature=24, continuation_steps=8, tolerance=1e-10,
           maximum_phase=.8):
    for value in (pump_frequency, signal_frequency, maximum_phase, tolerance):
        positive(value, 'frequencies/solver controls')
    positive(pump_power, 'pump_power', True)
    for value, name in [(harmonics, 'harmonics'), (sidebands, 'sidebands'), (continuation_steps, 'continuation_steps')]:
        if not isinstance(value, int) or value < 1:
            raise ValueError(f'{name} must be a positive integer')
    if not isinstance(quadrature, int) or quadrature < 8:
        raise ValueError('quadrature must be an integer >= 8')
    fi = pump_frequency-signal_frequency
    if min(signal_frequency, fi) <= 0 or abs(signal_frequency-fi)<pump_frequency*1e-10:
        raise ValueError('distinct positive 3WM signal and idler required')
    offset = signal_frequency-pump_frequency/2
    if sidebands*abs(offset) >= pump_frequency/2:
        raise ValueError('retained offset bands must stay within half the pump frequency')
    indices = [(2, 0), (1, 1), (1, -1)]
    indices += [(h, m) for h in range(0, 2*harmonics+1) for m in range(-sidebands, sidebands+1)
                if (h+m)%2==0 and h*pump_frequency/2+m*offset>0 and (h, m) not in indices]
    indices = tuple(indices)
    f = np.array([h*pump_frequency/2+m*(signal_frequency-pump_frequency/2) for h, m in indices])
    if np.any(f<=0) or np.min(np.diff(np.sort(f)))<pump_frequency*1e-10:
        raise ValueError('retained 3WM frequencies overlap or are nonpositive')
    _basis(indices, quadrature)
    c = chain.linear_chain
    potential = _Potential(chain.element)
    matrices = [_linear_system(c, frequency) for frequency in f]
    source = np.zeros((len(f), chain.cells+1), complex)
    source[0, 0] = 1j*np.sqrt(8*c.port_impedance*pump_power)/(c.port_impedance*c.critical_current)
    fields = np.zeros_like(source)
    for fraction in np.linspace(0, 1, continuation_steps+1)[1:]:
        fields, residual, _ = _newton(fields, matrices, source*fraction,
            tolerance=tolerance, nonlinearity=potential, indices=indices, quadrature=quadrature)
    _phase_check(fields, maximum_phase)
    return c, f, indices, matrices, source, fields, residual, potential


def _phase_check(fields, maximum):
    phase = float(np.max(np.sum(abs(np.diff(fields, axis=1)), axis=0)))
    if phase > maximum:
        raise ValueError(f'AC junction/element phase {phase:.3g} exceeds guard {maximum}')
    return phase


@dataclass
class ThreeWaveResult:
    frequency: np.ndarray
    mode_indices: tuple
    node_phase: np.ndarray
    dc_phase: np.ndarray
    gain: float
    idler_conversion: float
    pump_transmission: float
    transmitted_power: np.ndarray
    reflected_power: np.ndarray
    incident_power: np.ndarray
    maximum_phase: float
    residual: float
    depleted: bool


def three_wave_response(chain, pump_frequency, pump_power, signal_frequency,
                        signal_power=1e-20, *, depleted=False, **controls):
    positive(signal_power, 'signal_power')
    c, f, indices, matrices, source, fields, residual, potential = _state(
        chain, pump_frequency, pump_power, signal_frequency, **controls)
    q = controls.get('quadrature', 24)
    source[1, 0] = 1j*np.sqrt(8*c.port_impedance*signal_power)/(c.port_impedance*c.critical_current)
    if depleted:
        fields, residual, _ = _newton(fields, matrices, source,
            tolerance=controls.get('tolerance', 1e-10), nonlinearity=potential, indices=indices, quadrature=q)
    else:
        probe = np.zeros_like(source); probe[1] = source[1]
        fields += _complex_vector(spsolve(_jacobian(fields, matrices, potential, indices, q),
                                  _real_vector(probe)), chain.cells+1, len(f))
    phase = _phase_check(fields, controls.get('maximum_phase', .8))
    dc, _, _ = potential.dc(np.real(_basis(indices, q).T@(fields[:, :-1]-fields[:, 1:])))
    v = -2j*np.pi*f[:, None]*PHIBAR*fields
    incident = source[:, 0]*c.port_impedance*c.critical_current/2
    gain = float(abs(v[1, -1]/incident[1])**2)
    conversion = float(abs(v[2, -1]/incident[1])**2*signal_frequency/(pump_frequency-signal_frequency))
    pump = float(abs(v[0, -1]/incident[0])**2) if pump_power else np.nan
    return ThreeWaveResult(f, indices, fields, dc, gain, conversion, pump,
        abs(v[:, -1])**2/(2*c.port_impedance), abs(v[:, 0]-incident)**2/(2*c.port_impedance),
        abs(incident)**2/(2*c.port_impedance), phase, residual, depleted)


def quantum_noise(chain, pump_frequency, pump_power, signal_frequency, *,
                   port_temperature=0., junction_temperature=0., ground_temperature=0.,
                   resonator_temperature=0., commutator_tolerance=1e-6, **controls):
    """Pump-only quantum Floquet sidebands. Signal/idler are the first two modes.

    Excludes pump-frequency/DC sectors; no strong-probe quantum linearization.
    SNAIL temperature is common to all physical JJs in each cell.
    """
    from .jtwpa_noise import _noise_from_state
    c, f, indices, matrices, _, fields, residual, potential = _state(
        chain, pump_frequency, pump_power, signal_frequency, **controls)
    selected = [j for j, (harmonic, offset) in enumerate(indices) if harmonic % 2 and abs(offset)==1]
    jacobian = _jacobian(fields, matrices, potential, indices, controls.get('quadrature', 24))
    nodes, size = chain.cells+1, fields.size
    cols = np.concatenate([np.arange(j*nodes, (j+1)*nodes) for j in selected])
    both = np.r_[cols, size+cols]
    restricted = jacobian[both][:, both].tocsc()
    return _noise_from_state(c, f[selected], tuple(indices[j] for j in selected), fields[selected],
        restricted, residual, port_temperature=port_temperature,
        junction_temperature=junction_temperature, ground_temperature=ground_temperature,
        resonator_temperature=resonator_temperature, commutator_tolerance=commutator_tolerance,
        signal_mode=0)


def chain_netlist(chain, tones, transient=Transient(step=1e-12, stop=400e-9,
                   start=300e-9, output_step=10e-12), *, bias_ramp=100e-9,
                   snail_loop_inductance=1e-13, coupling=.2, coil_inductance=1e-9,
                   bias_shape='smooth'):
    """Physical full-sine RCSJ loops with inductors and DC flux coils.

    The SNAIL reduced model neglects geometric inductance; its netlist requires
    Lgeom/Lj_large <= 1e-3. Benchmark convergence toward that limit explicitly.
    """
    e = chain.element
    positive(bias_ramp, 'bias_ramp')
    positive(coil_inductance, 'coil_inductance')
    if not np.isfinite(coupling) or not 0<coupling<1:
        raise ValueError('coupling must lie between zero and one')
    if bias_shape not in ('smooth', 'linear'):
        raise ValueError('bias_shape must be smooth or linear')
    if isinstance(e, SNAIL):
        positive(snail_loop_inductance, 'snail_loop_inductance')
        if snail_loop_inductance*e.critical_current/PHIBAR > 1e-3:
            raise ValueError('SNAIL geometric inductance exceeds the supported reduced-model limit')
    tones = tuple(tones)
    if not tones:
        raise ValueError('at least one tone required')
    lines = ['* TWPA Sim: flux-biased 3WM physical RCSJ chain']
    lines += voltage_sources(tones, chain.port_impedance, prefix='VIN')
    lines += [f'RSOURCE SRC0 N0 {chain.port_impedance:.17g}',
              f'RLOAD N{chain.cells} 0 {chain.port_impedance:.17g}']
    def model(name, scale):
        return (f'.model {name} jj(RTYPE=0 IC={e.critical_current*scale:.17g} '
                f'C={e.junction_capacitance*scale:.17g} RN={e.junction_resistance/scale:.17g} '
                f'R0={e.junction_resistance/scale:.17g})')
    lines.append(model('JJL', 1.))
    if isinstance(e, SNAIL):
        lines.append(model('JJS', e.asymmetry))
    for cell in range(chain.cells):
        left, right = f'N{cell}', f'N{cell+1}'
        if isinstance(e, RFSQUID):
            lines.append(f'BJ{cell} {left} {right} JJL')
            loop, loop_start = e.loop_inductance, left
        else:
            lines.append(f'BJS{cell} {left} {right} JJS')
            previous = left
            for junction in range(e.junctions):
                target = f'Q{cell}_{junction}'
                lines.append(f'BJL{cell}_{junction} {previous} {target} JJL')
                previous = target
            loop, loop_start = snail_loop_inductance, previous
        mutual = coupling*np.sqrt(loop*coil_inductance)
        next_coil = f'COIL{cell+1}' if cell+1<chain.cells else '0'
        lines += [f'LLOOP{cell} {loop_start} {right} {loop:.17g}',
                  f'LCOIL{cell} COIL{cell} {next_coil} {coil_inductance:.17g}',
                  f'KBIAS{cell} LLOOP{cell} LCOIL{cell} {coupling:.17g}']
    # Series coils carry one prescribed current; their AC current is zero,
    # so the ideal bias network does not mediate inter-cell AC coupling.
    grid = np.linspace(0, 1, 257) if bias_shape=='smooth' else np.array([0., 1.])
    profile = 10*grid**3-15*grid**4+6*grid**5 if bias_shape=='smooth' else grid
    bias = e.flux_bias*PHI0/mutual
    points = ' '.join(f'{t*bias_ramp:.17g} {v*bias:.17g}' for t, v in zip(grid, profile))
    lines.append(f'IBIAS0 0 COIL0 PWL({points})')
    for node in range(chain.cells+1):
        weight = .5 if node in (0, chain.cells) else 1.
        lines.append(f'CG{node} N{node} 0 {chain.ground_capacitance*weight:.17g}')
        if chain.ground_resistance is not None:
            lines.append(f'RG{node} N{node} 0 {chain.ground_resistance/weight:.17g}')
        if chain.resonator is not None:
            r = chain.resonator
            lines += [f'CC{node} N{node} R{node} {r.coupling_capacitance*weight:.17g}',
                      f'CR{node} R{node} 0 {r.capacitance*weight:.17g}',
                      f'LR{node} R{node} 0 {r.inductance/weight:.17g}',
                      f'RR{node} R{node} 0 {r.resistance/weight:.17g}']
    middle = max(1, chain.cells//2)
    probes = sorted(set([0, middle-1, middle, min(chain.cells, middle+1), chain.cells]))
    trace = ' '.join(f'V(N{j})' for j in probes)
    trace += f' P(BJ{middle})' if isinstance(e, RFSQUID) else f' P(BJS{middle}) P(BJL{middle}_0)'
    lines += [f'.tran {transient.step:.17g} {transient.stop:.17g} {transient.start:.17g} {transient.output_step:.17g}',
              f'.print V(SRC0) {trace}', '.end']
    return '\n'.join(lines)+'\n'


def example_chain(kind='snail', cells=128, *, junction_resistance=1e9, ground_resistance=None,
                   dispersion_engineered=False):
    """Synthetic matched-circuit benchmark, not a fabricated-device preset."""
    if kind == 'rf-squid':
        trial = RFSQUID(1e-7, .8*PHIBAR/1e-7, 1e-12, junction_resistance, .32)
        ic = PHIBAR/(2.5e-9*coefficients(trial)[0])
        element = RFSQUID(ic, .8*PHIBAR/ic, 1e-12, junction_resistance, .32)
    elif kind == 'snail':
        trial = SNAIL(1e-7, .3, 3, 1e-12, junction_resistance, .35)
        ic = PHIBAR/(2.5e-9*coefficients(trial)[0])
        element = SNAIL(ic, .3, 3, 1e-12/(.3+1/3), junction_resistance, .35)
    else:
        raise ValueError('kind must be rf-squid or snail')
    chain = ThreeWaveChain(element, 1e-12, cells, ground_resistance=ground_resistance)
    from dataclasses import replace
    if dispersion_engineered:
        # Linear k_p = k_s+k_i at 1/.45/.55 GHz, with the unwanted sum bands
        # above the JJ-ladder cutoff. Values are not fitted to JoSIM gain.
        return replace(chain, element=replace(element, junction_capacitance=5*element.junction_capacitance),
            port_impedance=50., resonator=ResonantLoad(.2e-12, 1e-12,
            1/((2*np.pi*1.0584289309558132e9)**2*1e-12)))
    return replace(chain, port_impedance=float(bloch(chain, [1e9]).impedance[0].real))
