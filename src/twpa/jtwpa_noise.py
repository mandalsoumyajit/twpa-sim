"""Linearized quantum scattering of a pump-driven finite Josephson chain.

Photon-flux operators obey b = U a + V a^dagger + Ub c + Vb c^dagger.
Independent equilibrium resistor baths implement fluctuation-dissipation.
Spectral symmetrized noise is in quanta, with vacuum equal to one half.
"""
from dataclasses import dataclass
import numpy as np
from scipy.constants import h, k as kb
from scipy.sparse.linalg import splu
from .jpa_circuit import PHIBAR
from .jtwpa import _pump_solution, _jacobian, _newton
from .models import positive


def _temperatures(value, shape, name):
    try:
        result = np.broadcast_to(np.asarray(value, float), shape).copy()
    except ValueError as exc:
        raise ValueError(f'{name} must broadcast to {shape}') from exc
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise ValueError(f'{name} must be finite and nonnegative')
    return result


def _occupation(frequency, temperature):
    f, t = np.broadcast_arrays(frequency, temperature)
    ratio = np.divide(h*f, kb*t, out=np.full_like(f, np.inf), where=t>0)
    with np.errstate(over='ignore'):
        return 1/np.expm1(ratio)


def _covariances(u, v, variance):
    return ((u*variance)@u.conj().T+(v*variance)@v.conj().T,
            (u*variance)@v.T+(v*variance)@u.T)


def _quadratures(normal, anomalous):
    return np.block([[(normal+anomalous).real, (anomalous-normal).imag],
                     [(anomalous+normal).imag, (normal-anomalous).real]])


@dataclass
class QuantumNoiseResult:
    frequency: np.ndarray
    mode_indices: tuple
    node_phase: np.ndarray  # Coherent operating point, including probe if depleted.
    annihilation: np.ndarray  # Output vs incoming port annihilation operators.
    creation: np.ndarray
    bath_annihilation: np.ndarray
    bath_creation: np.ndarray
    bath_labels: tuple  # (kind, frequency index, site index) per bath column.
    port_occupation: np.ndarray
    bath_occupation: np.ndarray
    covariance: np.ndarray  # <{b_i,b_j^dagger}>/2, vacuum I/2.
    anomalous_covariance: np.ndarray  # <{b_i,b_j}>/2.
    quadrature_covariance: np.ndarray  # [all x, all p], [x,p]=i.
    covariance_contributions: dict  # Ports and each resistor family.
    gain: float  # Right signal output vs left signal input, |U|^2.
    output_noise: float
    added_noise: float  # Input referred, subtract incoming signal n+1/2.
    quantum_limit: float  # |1-1/G|/2, phase-preserving vacuum bound.
    excess_noise: float
    phase_preserving: bool
    commutator_error: float
    anomalous_commutator_error: float
    minimum_uncertainty_eigenvalue: float
    residual: float

    @property
    def gain_db(self):
        return 10*np.log10(self.gain)

    @property
    def output_photons(self):
        """Normally ordered spectral occupation, excluding coherent pump tones."""
        return np.diag(self.covariance).real-.5


def four_wave_noise(chain, pump_frequency, pump_power, signal_frequency, *,
                    signal_power=0., depleted=False,
                    port_temperature=0., junction_temperature=0.,
                    ground_temperature=0., resonator_temperature=0.,
                    continuation_steps=8, tolerance=1e-10, maximum_phase=1.2,
                    nonlinearity='sine', harmonics=None, sidebands=1,
                    quadrature=16, commutator_tolerance=1e-6):
    """Gaussian quantum noise about a classical operating point.

    Port channels are frequency-major (left, right). Internal bath columns
    are identified by bath_labels; each physical resistor has an independent
    input at each retained frequency. Temperatures are kelvin: ports broadcast
    to (modes,2), junctions to (cells,), other baths to (cells+1,).

    By default the operating point is pump only. With depleted=True and a
    positive signal_power, first solve the coherent multitone state. Scalar
    added noise is NaN when that state's signal response is phase sensitive.
    This is not quantum pump depletion or a dynamic stability test.
    A noncanonical transfer raises rather than being repaired.
    """
    positive(commutator_tolerance, 'commutator_tolerance')
    positive(signal_power, 'signal_power', True)
    if bool(depleted) != bool(signal_power):
        raise ValueError('depleted=True requires positive signal_power, otherwise use zero')
    f, indices, matrices, source, fields, residual = _pump_solution(
        chain, pump_frequency, pump_power, signal_frequency,
        continuation_steps=continuation_steps, tolerance=tolerance,
        maximum_phase=maximum_phase, nonlinearity=nonlinearity,
        harmonics=harmonics, sidebands=sidebands, quadrature=quadrature)
    if depleted:
        source[1, 0] = 1j*np.sqrt(8*chain.port_impedance*signal_power)/(chain.port_impedance*chain.critical_current)
        fields, residual, _ = _newton(fields, matrices, source,
            tolerance=tolerance, nonlinearity=nonlinearity, indices=indices, quadrature=quadrature)
        phase = float(np.max(np.sum(abs(fields[:, :-1]-fields[:, 1:]), axis=0)))
        if phase > maximum_phase:
            raise ValueError(f'Junction phase {phase:.3g} exceeds the retained-mode guard {maximum_phase}')
    modes, nodes = fields.shape
    ports = 2*modes
    size = modes*nodes
    tp = _temperatures(port_temperature, (modes, 2), 'port_temperature')
    tj = _temperatures(junction_temperature, (chain.cells,), 'junction_temperature')
    tg = _temperatures(ground_temperature, (nodes,), 'ground_temperature')
    tr = _temperatures(resonator_temperature, (nodes,), 'resonator_temperature')

    # Selected-output adjoint solves avoid a dense inverse or one solve for
    # each of the O(cells*modes) internal bath channels.
    select = np.zeros((2*ports, 2*size))
    scale = 2*np.pi*f*PHIBAR/np.sqrt(2*chain.port_impedance*h*f)
    for mode in range(modes):
        for port, node in enumerate((0, nodes-1)):
            row, col = 2*mode+port, mode*nodes+node
            select[row, size+col] = scale[mode]
            select[ports+row, col] = -scale[mode]
    factor = splu(_jacobian(fields, matrices, nonlinearity, indices, quadrature))
    transfer = factor.solve(select.T, trans='T').T
    output = transfer[:ports]+1j*transfer[ports:]

    def couple(mode, sites, conductance, other=None):
        sites = np.atleast_1d(sites)
        cols = mode*nodes+sites
        real, imag = output[:, cols], output[:, size+cols]
        if other is not None:
            real = real-output[:, mode*nodes+other]
            imag = imag-output[:, size+mode*nodes+other]
        source = np.sqrt(8*h*f[mode]*conductance)/chain.critical_current
        return (real-1j*imag)*source/2, (real+1j*imag)*source/2

    u, v = [], []
    for mode in range(modes):
        a, b = couple(mode, [0, nodes-1], 1/chain.port_impedance)
        u.append(a); v.append(b)
    u, v = np.concatenate(u, axis=1)-np.eye(ports), np.concatenate(v, axis=1)

    weights = np.ones(nodes); weights[[0, -1]] = .5
    groups = [('junction', np.arange(chain.cells), np.arange(1, nodes),
               np.full((modes, chain.cells), 1/chain.junction_resistance), tj)]
    if chain.ground_resistance is not None:
        groups.append(('ground', np.arange(nodes), None,
                       np.broadcast_to(weights/chain.ground_resistance, (modes, nodes)), tg))
    if chain.resonator is not None:
        r = chain.resonator
        w = 2*np.pi*f
        yc = -1j*w*r.coupling_capacitance
        yp = 1/r.resistance-1j*w*r.capacitance+1/(-1j*w*r.inductance)
        # Eliminate the resonator node: G_eff = |Yc/(Yc+Yparallel)|^2/R.
        # Its phase can be absorbed into each independent thermal bath input.
        conductance = abs(yc/(yc+yp))**2/r.resistance
        groups.append(('resonator', np.arange(nodes), None,
                       conductance[:, None]*weights, tr))
    ub, vb, occupations, labels, slices = [], [], [], [], {}
    offset = 0
    for kind, sites, other, conductance, temperature in groups:
        for mode in range(modes):
            a, b = couple(mode, sites, conductance[mode], other)
            ub.append(a); vb.append(b)
            occupations.extend(_occupation(f[mode], temperature))
            labels.extend((kind, mode, int(site)) for site in sites)
        count = modes*len(sites)
        slices[kind] = slice(offset, offset+count)
        offset += count
    ub, vb = np.concatenate(ub, axis=1), np.concatenate(vb, axis=1)
    nb = np.array(occupations)
    np_ports = _occupation(f[:, None], tp).ravel()
    normal, anomalous = _covariances(u, v, np_ports+.5)
    contributions = {'ports': normal.copy()}
    for kind, columns in slices.items():
        cn, cm = _covariances(ub[:, columns], vb[:, columns], nb[columns]+.5)
        normal += cn; anomalous += cm
        contributions[kind] = cn
    canonical = u@u.conj().T-v@v.conj().T+ub@ub.conj().T-vb@vb.conj().T
    pair = u@v.T-v@u.T+ub@vb.T-vb@ub.T
    error = float(np.max(abs(canonical-np.eye(ports))))
    pair_error = float(np.max(abs(pair)))
    if max(error, pair_error) > commutator_tolerance:
        raise RuntimeError(f'Quantum commutator defect {max(error, pair_error):.3g}; '
                           'refine the pump solve/Fourier quadrature')
    covariance = _quadratures(normal, anomalous)
    omega = np.block([[np.zeros((ports, ports)), np.eye(ports)],
                      [-np.eye(ports), np.zeros((ports, ports))]])
    minimum = float(np.linalg.eigvalsh(covariance+.5j*omega)[0])
    incoming, outgoing = 2, 3  # Primary signal, left -> right.
    gain = float(abs(u[outgoing, incoming])**2)
    phase_preserving = bool(abs(v[outgoing, incoming]) < 1e-8*max(1., np.sqrt(gain)))
    noise = float(normal[outgoing, outgoing].real)
    added = float(noise/gain-(np_ports[incoming]+.5)) if gain>0 and phase_preserving else np.nan
    limit = float(abs(1-1/gain)/2) if gain>0 and phase_preserving else np.nan
    return QuantumNoiseResult(f, indices, fields, u, v, ub, vb, tuple(labels),
                np_ports, nb, normal, anomalous, covariance, contributions,
                gain, noise, added, limit, added-limit, phase_preserving,
                error, pair_error, minimum, residual)
