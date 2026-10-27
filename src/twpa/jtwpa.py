"""Finite Josephson-chain 4WM harmonic-balance models.

Peak flux phasors use exp(-i omega t). Exact-sine Fourier projection or a
quartic-potential approximation includes 2*pump <-> signal+idler mixing.
Optional odd-harmonic bands and sidebands expose truncation errors. Finite
circuit KCL includes boundaries, loss and optional resonant phase matching.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.sparse import diags, bmat
from scipy.sparse.linalg import spsolve
from .jpa_circuit import PHIBAR, Tone, Transient, voltage_sources
from .models import positive
from .linear import _power


@dataclass(frozen=True)
class ResonantLoad:
    coupling_capacitance: float
    capacitance: float
    inductance: float
    resistance: float = 1e9  # Parallel resonator loss.

    def __post_init__(self):
        for key, value in vars(self).items():
            positive(value, key)

    def admittance(self, frequency):
        w = 2*np.pi*np.asarray(frequency, float)
        result = np.zeros_like(w, complex)
        nz = w != 0
        yp = 1/self.resistance-1j*w[nz]*self.capacitance+1/(-1j*w[nz]*self.inductance)
        result[nz] = 1/(1/(-1j*w[nz]*self.coupling_capacitance)+1/yp)
        return result


@dataclass(frozen=True)
class Chain:
    critical_current: float
    ground_capacitance: float
    junction_capacitance: float
    cells: int
    port_impedance: float = 50.
    junction_resistance: float = 1e9
    ground_resistance: float | None = None
    cell_length: float = 10e-6
    resonator: ResonantLoad | None = None

    def __post_init__(self):
        for key in ('critical_current', 'ground_capacitance', 'junction_capacitance',
                    'port_impedance', 'junction_resistance', 'cell_length'):
            positive(getattr(self, key), key)
        if not isinstance(self.cells, int) or self.cells < 1:
            raise ValueError('cells must be a positive integer')
        if self.ground_resistance is not None:
            positive(self.ground_resistance, 'ground_resistance')

    @property
    def inductance(self):
        return PHIBAR/self.critical_current

    @property
    def plasma_frequency(self):
        return 1/(2*np.pi*np.sqrt(self.inductance*self.junction_capacitance))

    def admittances(self, frequency):
        f = np.atleast_1d(np.asarray(frequency, float))
        if np.any(~np.isfinite(f)) or np.any(f < 0):
            raise ValueError('frequency must be finite and nonnegative')
        w = 2*np.pi*f
        series = np.zeros_like(w, complex)
        nz = w != 0
        series[nz] = 1/(1/self.junction_resistance-1j*w[nz]*self.junction_capacitance+
                        1/(-1j*w[nz]*self.inductance))
        shunt = -1j*w*self.ground_capacitance
        if self.ground_resistance is not None:
            shunt += 1/self.ground_resistance
        if self.resonator is not None:
            shunt += self.resonator.admittance(f)
        return series, shunt

    def cell_matrix(self, frequency, *, normalized=False):
        z, y = self.admittances(frequency)
        t = np.zeros((len(z), 2, 2), complex)
        t[:, 0, 0] = t[:, 1, 1] = 1+z*y/2
        t[:, 0, 1] = z
        t[:, 1, 0] = y*(1+z*y/4)
        if normalized:
            t[:, 0, 1] /= self.port_impedance
            t[:, 1, 0] *= self.port_impedance
        return t


@dataclass
class BlochResult:
    frequency: np.ndarray
    k: np.ndarray  # Principal-zone radians/cell; use finite scattering at band edges.
    attenuation: np.ndarray
    impedance: np.ndarray
    passband: np.ndarray


def bloch(chain, frequency):
    f = np.atleast_1d(np.asarray(frequency, float))
    t = chain.cell_matrix(f)
    halftrace = (t[:, 0, 0]+t[:, 1, 1])/2
    theta = np.arccos(halftrace.astype(complex))
    impedance = np.zeros_like(theta)
    nz = f != 0
    impedance[nz] = np.sqrt(t[nz, 0, 1]/t[nz, 1, 0])
    impedance[impedance.real < 0] *= -1
    c0 = chain.ground_capacitance
    if chain.resonator is not None:
        c0 += chain.resonator.coupling_capacitance
    impedance[~nz] = np.sqrt(chain.inductance/c0) if chain.ground_resistance is None else 0.
    w = 2*np.pi*f
    denominator = 1-w*w*chain.inductance*chain.junction_capacitance
    ceff = np.full_like(f, chain.ground_capacitance)
    poles = abs(denominator) < 1e-12
    if chain.resonator is not None:
        r = chain.resonator
        d = 1-w*w*r.inductance*(r.capacitance+r.coupling_capacitance)
        poles |= abs(d)<1e-12
        ceff += np.divide(r.coupling_capacitance*(1-w*w*r.inductance*r.capacitance),
                          d, out=np.full_like(f, np.inf), where=abs(d)>=1e-12)
    ideal_trace = 1-np.divide(w*w*chain.inductance*ceff, 2*denominator,
                              out=np.full_like(f, np.inf), where=abs(denominator)>=1e-12)
    return BlochResult(f, theta.real, abs(theta.imag), impedance,
                       (abs(ideal_trace)<=1+1e-12)&~poles)


def linear_response(chain, frequency):
    """Exact finite-network scattering, including half shunts at both ends."""
    t = chain.cell_matrix(frequency, normalized=True)
    t, scale = _power(t, np.zeros(len(t)), chain.cells)
    denominator = t.sum(axis=(1, 2))
    return dict(s21=2*np.exp(-scale)/denominator,
                s11=(t[:, 0, 0]+t[:, 0, 1]-t[:, 1, 0]-t[:, 1, 1])/denominator,
                s21_db=20/np.log(10)*(np.log(2)-scale-np.log(abs(denominator))))


def _cube(phasors):
    """Peak phasors of (Re[P ep+S es+I ei])^3 at [p,s,i]."""
    p, s, i = phasors
    return np.array([
        .75*abs(p)**2*p+1.5*(abs(s)**2+abs(i)**2)*p+1.5*s*i*np.conj(p),
        .75*abs(s)**2*s+1.5*(abs(p)**2+abs(i)**2)*s+.75*p*p*np.conj(i),
        .75*abs(i)**2*i+1.5*(abs(p)**2+abs(s)**2)*i+.75*p*p*np.conj(s)])


def _linear_system(chain, frequency):
    n = chain.cells+1
    w = 2*np.pi*frequency
    l = chain.inductance
    branch = 1-w*w*l*chain.junction_capacitance-1j*w*l/chain.junction_resistance
    weights = np.ones(n); weights[[0, -1]] = .5
    _, shunt = chain.admittances([frequency])
    diagonal = weights*(-1j*w*l*shunt[0])
    diagonal[0] += branch-1j*w*l/chain.port_impedance
    diagonal[-1] += branch-1j*w*l/chain.port_impedance
    diagonal[1:-1] += 2*branch
    return diags([-branch*np.ones(n-1), diagonal, -branch*np.ones(n-1)], [-1, 0, 1], format='csc')


def _real_vector(phasors):
    return np.r_[phasors.real.ravel(), phasors.imag.ravel()]


def _complex_vector(vector, nodes, modes=3):
    return (vector[:modes*nodes]+1j*vector[modes*nodes:]).reshape(modes, nodes)


def _incidence(values):
    return np.pad(values, ((0, 0), (0, 1)))-np.pad(values, ((0, 0), (1, 0)))


@lru_cache(maxsize=16)
def _basis(indices, quadrature):
    if quadrature < 2*max(max(abs(h), abs(m)) for h, m in indices)+4:
        raise ValueError('spectral quadrature is too small for retained modes')
    a, b = np.meshgrid(np.arange(quadrature)*2*np.pi/quadrature,
                       np.arange(quadrature)*2*np.pi/quadrature, indexing='ij')
    return np.array([np.exp(-1j*(h*a+m*b)).ravel() for h, m in indices])


def _nonlinear(branch, nonlinearity='cubic', indices=None, quadrature=16, derivative=False):
    modes = len(branch)
    if nonlinearity == 'cubic':
        if modes != 3:
            raise ValueError('cubic closed form requires the three fundamental modes')
        current = _cube(branch)/6
        if not derivative:
            return current
        derivatives = []
        step = 1e-6
        for col in range(2*modes):
            perturb = np.zeros_like(branch)
            perturb[col % modes] = step*(1 if col < modes else 1j)
            derivatives.append((_cube(branch+perturb)-_cube(branch-perturb))/(12*step))
        return derivatives
    basis = _basis(tuple(indices), quadrature)
    theta = np.real(basis.T@branch)
    projection = 2*basis.conj()/basis.shape[1]
    if not derivative:
        return branch-projection@np.sin(theta)
    cosine = np.cos(theta)
    derivatives = []
    for col in range(2*modes):
        local = basis[col % modes].real if col < modes else -basis[col % modes].imag
        d = -projection@(cosine*local[:, None])
        d[col % modes] += 1 if col < modes else 1j
        derivatives.append(d)
    return derivatives


def _residual(fields, matrices, source, nonlinearity='cubic', indices=None, quadrature=16):
    branch = fields[:, :-1]-fields[:, 1:]
    return np.array([a@v for a, v in zip(matrices, fields)])-_incidence(
        _nonlinear(branch, nonlinearity, indices, quadrature))-source


def _jacobian(fields, matrices, nonlinearity='cubic', indices=None, quadrature=16):
    n = fields.shape[1]
    modes = len(fields)
    branch = fields[:, :-1]-fields[:, 1:]
    # Cubic current is local to each branch. Vectorized local derivatives
    # assemble a sparse analytic-circuit Jacobian, without dense global differences.
    derivatives = _nonlinear(branch, nonlinearity, indices, quadrature, derivative=True)
    def laplace(values):
        diagonal = np.r_[values, 0]+np.r_[0, values]
        return diags([-values, diagonal, -values], [-1, 0, 1], shape=(n, n), format='csc')
    blocks = []
    for row in range(2*modes):
        sub = []
        for col in range(2*modes):
            component = derivatives[col][row % modes]
            nonlinear = laplace(component.real if row < modes else component.imag)
            if row % modes == col % modes:
                a = matrices[row % modes]
                linear = a.real if (row<modes)==(col<modes) else (-a.imag if row<modes else a.imag)
                sub.append(linear-nonlinear)
            else:
                sub.append(-nonlinear)
        blocks.append(sub)
    return bmat(blocks, format='csc')


def _newton(initial, matrices, source, tolerance=1e-10, iterations=35,
            nonlinearity='cubic', indices=None, quadrature=16):
    fields = initial.copy()
    denominator = max(float(np.max(abs(source))), 1e-12)
    for iteration in range(iterations):
        residual = _residual(fields, matrices, source, nonlinearity, indices, quadrature)
        norm = float(np.max(abs(residual)))/denominator
        if norm < tolerance:
            return fields, norm, iteration
        step = _complex_vector(spsolve(_jacobian(fields, matrices, nonlinearity, indices, quadrature),
                                       -_real_vector(residual)), fields.shape[1], len(fields))
        damping = 1.
        for _ in range(14):
            candidate = fields+damping*step
            candidate_norm = float(np.max(abs(_residual(candidate, matrices, source, nonlinearity, indices, quadrature))))/denominator
            if candidate_norm < norm:
                fields = candidate
                break
            damping *= .5
        else:
            raise RuntimeError('4WM Newton line search failed; reduce drive or refine continuation')
    raise RuntimeError('4WM harmonic balance did not converge')


@dataclass
class FourWaveResult:
    frequency: np.ndarray  # First pump, signal, idler, then extra retained modes [Hz].
    node_phase: np.ndarray
    gain: float
    idler_conversion: float  # photon flux, not power.
    pump_transmission: float
    reflection: complex
    residual: float
    maximum_junction_phase: float  # Conservative sum of branch-tone magnitudes.
    depleted: bool
    incident_power: np.ndarray
    transmitted_power: np.ndarray
    reflected_power: np.ndarray
    mode_indices: tuple
    nonlinearity: str


def _pump_solution(chain, pump_frequency, pump_power, signal_frequency,
                       *, continuation_steps=8,
                       tolerance=1e-10, maximum_phase=1.2, nonlinearity='sine',
                       harmonics=None, sidebands=1, quadrature=16):
    """Shared pump-only operating point for classical and quantum linearization."""
    for value in (pump_frequency, signal_frequency, maximum_phase, tolerance):
        positive(value, 'frequencies/power/solver controls')
    positive(pump_power, 'pump_power', True)
    if not isinstance(continuation_steps, int) or continuation_steps < 1:
        raise ValueError('continuation_steps must be a positive integer')
    if nonlinearity not in ('cubic', 'sine'):
        raise ValueError('nonlinearity must be cubic or sine')
    if harmonics is None:
        harmonics = 1 if nonlinearity == 'cubic' else 3
    if not isinstance(harmonics, int) or harmonics < 1 or harmonics % 2 == 0:
        raise ValueError('harmonics must be a positive odd integer')
    if not isinstance(sidebands, int) or sidebands < 1:
        raise ValueError('sidebands must be a positive integer')
    if not isinstance(quadrature, int) or quadrature < 8:
        raise ValueError('quadrature must be an integer >= 8')
    if nonlinearity == 'cubic' and (harmonics != 1 or sidebands != 1):
        raise ValueError('cubic closed form retains three tones only')
    fi = 2*pump_frequency-signal_frequency
    if fi <= 0 or abs(fi-signal_frequency) < pump_frequency*1e-10:
        raise ValueError('distinct positive signal and 4WM idler frequencies required')
    indices = [(1, 0), (1, 1), (1, -1)]
    indices += [(h, m) for h in range(1, harmonics+1, 2) for m in range(-sidebands, sidebands+1)
                if (h, m) not in indices]
    indices = tuple(indices)
    frequencies = np.array([h*pump_frequency+m*(signal_frequency-pump_frequency) for h, m in indices])
    # Exclude accidental third-harmonic aliases in the three-tone approximation.
    if max(frequencies[:3])/min(frequencies[:3]) >= 2 or np.any(frequencies <= 0) or len(np.unique(frequencies)) != len(frequencies):
        raise ValueError('three-tone model requires nearby pump, signal and idler bands')
    matrices = [_linear_system(chain, f) for f in frequencies]
    source = np.zeros((len(indices), chain.cells+1), complex)
    source[0, 0] = 1j*np.sqrt(8*chain.port_impedance*pump_power)/(chain.port_impedance*chain.critical_current)
    pump_source = np.zeros_like(source); pump_source[0] = source[0]
    fields = np.zeros_like(source)
    for fraction in np.linspace(0, 1, continuation_steps+1)[1:]:
        fields, residual, _ = _newton(fields, matrices, pump_source*fraction,
                 tolerance=tolerance, nonlinearity=nonlinearity, indices=indices, quadrature=quadrature)
    phase = float(max(np.sum(abs(fields[:, :-1]-fields[:, 1:]), axis=0)))
    if phase > maximum_phase:
        raise ValueError(f'Junction phase {phase:.3g} exceeds the retained-mode guard {maximum_phase}')
    return frequencies, indices, matrices, source, fields, residual


def four_wave_response(chain, pump_frequency, pump_power, signal_frequency,
                       signal_power=1e-20, *, depleted=False, continuation_steps=8,
                       tolerance=1e-10, maximum_phase=1.2, nonlinearity='sine',
                       harmonics=None, sidebands=1, quadrature=16):
    """Multitone classical 4WM scattering at finite resistive ports.

    Undepleted mode linearizes about a self-consistent pump-only solution.
    Depleted mode solves all retained amplitudes together. Exact-sine Fourier
    projection is the default; 'cubic' retains a three-tone quartic potential.
    A converged algebraic solution does not establish dynamical stability.
    """
    positive(signal_power, 'signal_power')
    frequencies, indices, matrices, source, fields, residual = _pump_solution(
        chain, pump_frequency, pump_power, signal_frequency,
        continuation_steps=continuation_steps, tolerance=tolerance,
        maximum_phase=maximum_phase, nonlinearity=nonlinearity,
        harmonics=harmonics, sidebands=sidebands, quadrature=quadrature)
    fi = 2*pump_frequency-signal_frequency
    source[1, 0] = 1j*np.sqrt(8*chain.port_impedance*signal_power)/(chain.port_impedance*chain.critical_current)
    if depleted:
        fields, residual, _ = _newton(fields, matrices, source, tolerance=tolerance,
                  nonlinearity=nonlinearity, indices=indices, quadrature=quadrature)
    else:
        perturbation = np.zeros_like(source); perturbation[1] = source[1]
        correction = _complex_vector(spsolve(_jacobian(fields, matrices, nonlinearity, indices, quadrature),
                                             _real_vector(perturbation)), chain.cells+1, len(indices))
        fields += correction
        # Residual is for the pump solve; the probe is a first-order solution.
    phase = float(max(np.sum(abs(fields[:, :-1]-fields[:, 1:]), axis=0)))
    if phase > maximum_phase:
        raise ValueError(f'Junction phase {phase:.3g} exceeds the retained-mode guard {maximum_phase}')
    voltage = -2j*np.pi*frequencies[:, None]*PHIBAR*fields
    vs = source[:, 0]*chain.port_impedance*chain.critical_current
    incident = vs/2
    gain = float(abs(voltage[1, -1]/incident[1])**2)
    conversion = float(abs(voltage[2, -1]/incident[1])**2*signal_frequency/fi)
    pump_transmission = float(abs(voltage[0, -1]/incident[0])**2) if pump_power else np.nan
    reflection = (voltage[1, 0]-incident[1])/incident[1]
    incoming_power = abs(incident)**2/(2*chain.port_impedance)
    transmitted_power = abs(voltage[:, -1])**2/(2*chain.port_impedance)
    reflected_power = abs(voltage[:, 0]-incident)**2/(2*chain.port_impedance)
    return FourWaveResult(frequencies, fields, gain, conversion, pump_transmission,
                          complex(reflection), residual, phase, depleted,
                          incoming_power, transmitted_power, reflected_power, indices, nonlinearity)


def obrien_2014(cells=2000):
    """Published RPM circuit values, with Lj derived consistently from Ic.

    The paper quotes rounded L=100 pH and Ic=3.29 uA; phi_bar/Ic gives
    100.0337 pH. This is a circuit preset, not a claim of reproducing every
    published continuum-envelope curve with finite-port harmonic balance.
    """
    return Chain(3.29e-6, 39e-15, 329e-15, cells, cell_length=10e-6,
                 resonator=ResonantLoad(10e-15, 7.036e-12, 100e-12))


def quantum_noise(chain, pump_frequency, pump_power, signal_frequency, **kwargs):
    """Quantum small-signal noise; see twpa.jtwpa_noise.four_wave_noise."""
    from .jtwpa_noise import four_wave_noise
    return four_wave_noise(chain, pump_frequency, pump_power, signal_frequency, **kwargs)


def chain_netlist(chain, tones, transient=Transient(step=2e-12, stop=200e-9, start=100e-9, output_step=10e-12)):
    """RCSJ chain with symmetric half-shunts and physical matched-resistor ports."""
    tones = tuple(tones)
    if not tones:
        raise ValueError('provide at least one tone')
    lines = ['* TWPA Sim: finite Josephson junction chain, full sine CPR']
    lines += voltage_sources(tones, chain.port_impedance, prefix='VIN')
    lines += [f'RSOURCE SRC0 N0 {chain.port_impedance:.17g}',
              f'RLOAD N{chain.cells} 0 {chain.port_impedance:.17g}',
              f'.model JJMOD jj(RTYPE=0 IC={chain.critical_current:.17g} '
              f'C={chain.junction_capacitance:.17g} RN={chain.junction_resistance:.17g} '
              f'R0={chain.junction_resistance:.17g})']
    for i in range(chain.cells):
        lines.append(f'BJ{i} N{i} N{i+1} JJMOD')
    for i in range(chain.cells+1):
        factor = .5 if i in (0, chain.cells) else 1.
        lines.append(f'CG{i} N{i} 0 {chain.ground_capacitance*factor:.17g}')
        if chain.ground_resistance is not None:
            lines.append(f'RG{i} N{i} 0 {chain.ground_resistance/factor:.17g}')
        if chain.resonator is not None:
            r = chain.resonator
            lines += [f'CC{i} N{i} Q{i} {r.coupling_capacitance*factor:.17g}',
                      f'CR{i} Q{i} 0 {r.capacitance*factor:.17g}',
                      f'LR{i} Q{i} 0 {r.inductance/factor:.17g}',
                      f'RR{i} Q{i} 0 {r.resistance/factor:.17g}']
    middle = max(1, chain.cells//2)
    probes = sorted(set([0, middle-1, middle, min(chain.cells, middle+1), chain.cells]))
    trace = ' '.join(f'V(N{i})' for i in probes)
    lines += [f'.tran {transient.step:.17g} {transient.stop:.17g} {transient.start:.17g} {transient.output_step:.17g}',
              f'.print V(SRC0) {trace}', '.end']
    return '\n'.join(lines)+'\n'
