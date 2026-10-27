"""Single-mode resonant Josephson parametric amplifiers.

Convention: exp(-i omega t), a_out = sqrt(kappa_ext)*a - a_in.
Frequencies are Hz; kappa, Kerr K and pair drive are angular rates (rad/s).
|a|^2 is cavity occupation; |a_in|^2 is incident photons/second.
"""
from dataclasses import dataclass
import numpy as np
from scipy.constants import h, k as kb
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from .models import positive


@dataclass(frozen=True)
class Resonator:
    frequency: float
    kappa_external: float
    kappa_internal: float = 0.
    kerr: float = 0.
    internal_temperature: float = 0.

    def __post_init__(self):
        positive(self.frequency, 'frequency')
        positive(self.kappa_external, 'kappa_external')
        positive(self.kappa_internal, 'kappa_internal', True)
        positive(self.internal_temperature, 'internal_temperature', True)
        if not np.isfinite(self.kerr):
            raise ValueError('kerr must be finite')

    @property
    def kappa(self):
        return self.kappa_external+self.kappa_internal


def _drift(resonator, detuning, pair_drive):
    return np.array([[-resonator.kappa/2+1j*detuning, pair_drive],
                     [np.conj(pair_drive), -resonator.kappa/2-1j*detuning]])


def _stability(resonator, matrix):
    eigenvalues = np.linalg.eigvals(matrix)
    margin = -float(max(eigenvalues.real))
    tol = resonator.kappa*1e-9
    label = 'stable' if margin > tol else ('unstable' if margin < -tol else 'marginal')
    return eigenvalues, margin, label


@dataclass(frozen=True)
class PumpPoint:
    frequency: float
    power: float
    phase: float
    amplitude: complex
    photons: float
    eigenvalues: np.ndarray
    stability: str
    stability_margin: float
    resonator: Resonator


def pump_points(resonator, frequency, power, phase=0.):
    """All nonnegative Kerr pump roots, sorted by occupation, with stability.

    Bracket the dimensionless cubic between its turning points, retaining
    repeated roots at folds. No branch is silently selected in bistability.
    """
    positive(frequency, 'pump frequency')
    positive(power, 'pump power', True)
    if not np.isfinite(phase):
        raise ValueError('pump phase must be finite')
    r = resonator
    delta = 2*np.pi*(frequency-r.frequency)
    incoming = np.sqrt(power/(h*frequency))*np.exp(1j*phase)
    if r.kerr == 0 or power == 0:
        photons = [r.kappa_external*abs(incoming)**2/((r.kappa/2)**2+delta**2)]
    else:
        d, sign = delta/r.kappa, np.sign(r.kerr)
        drive = abs(r.kerr)*r.kappa_external*abs(incoming)**2/r.kappa**3
        def polynomial(x):
            return x*(.25+(d-sign*x)**2)-drive
        turns = []
        disc = d*d-.75
        if disc >= 0:
            turns = [x for x in ((2*d*sign-np.sqrt(disc))/3,
                                (2*d*sign+np.sqrt(disc))/3) if x > 0]
        upper = max(1., 2*abs(d)+2*np.cbrt(drive)+1)
        knots = sorted(set([0., *turns, upper]))
        roots = []
        tol = 64*np.finfo(float).eps*max(drive, np.finfo(float).tiny)
        for x in knots:
            if abs(polynomial(x)) <= tol:
                roots.append(x)
        for lo, hi in zip(knots[:-1], knots[1:]):
            flo, fhi = polynomial(lo), polynomial(hi)
            if abs(flo) <= tol or abs(fhi) <= tol:
                continue
            if np.signbit(flo) != np.signbit(fhi):
                root = brentq(polynomial, lo, hi, xtol=np.finfo(float).tiny, rtol=1e-14)
                if not any(abs(root-x) < 1e-9*max(1., abs(root)) for x in roots):
                    roots.append(root)
        photons = sorted(x*r.kappa/abs(r.kerr) for x in roots)
    points = []
    for n in photons:
        amplitude = np.sqrt(r.kappa_external)*incoming/(r.kappa/2-1j*(delta-r.kerr*n))
        matrix = _drift(r, delta-2*r.kerr*n, -1j*r.kerr*amplitude**2)
        eigenvalues, margin, label = _stability(r, matrix)
        points.append(PumpPoint(frequency, power, phase, complex(amplitude), float(n),
                                eigenvalues, label, margin, r))
    return points


def select_pump_point(points, branch=None, previous=None):
    """Choose a stable root explicitly, or follow nearest prior amplitude.

    'lower'/'upper' refer to occupation. Discrete continuation follows stable
    roots but does not model switching rates or dynamic hysteresis.
    """
    stable = [p for p in points if p.stability == 'stable']
    if not stable:
        raise ValueError('no strictly stable pump operating point')
    if branch not in (None, 'lower', 'upper'):
        raise ValueError('branch must be lower or upper')
    if branch is not None and previous is not None:
        raise ValueError('choose branch or continuation, not both')
    if previous is not None:
        return min(stable, key=lambda p: abs(p.amplitude-previous.amplitude))
    if branch == 'lower':
        return min(stable, key=lambda p: p.photons)
    if branch == 'upper':
        return max(stable, key=lambda p: p.photons)
    if len(stable) != 1:
        raise ValueError('multiple stable pump roots; choose a branch or previous point')
    return stable[0]


def _occupation(f, temperature):
    if temperature == 0:
        return np.zeros_like(f)
    return 1/np.expm1(np.minimum(h*f/(kb*temperature), 700))


@dataclass
class JPAResponse:
    signal_frequency: np.ndarray
    idler_frequency: np.ndarray
    frame_frequency: float
    scattering: np.ndarray  # [a_s, a_i^dagger] output vs external input
    bath_scattering: np.ndarray
    gain: np.ndarray
    conversion_gain: np.ndarray
    added_noise: np.ndarray  # NaN at degeneracy: use quadratures instead.
    degenerate: np.ndarray
    eigenvalues: np.ndarray
    resonator: Resonator

    @property
    def gain_db(self):
        return 10*np.log10(self.gain)


def _response(r, frame, detuning, pair_drive, signal_frequency):
    f = np.atleast_1d(np.asarray(signal_frequency, float))
    if f.ndim != 1 or not len(f) or np.any(~np.isfinite(f)):
        raise ValueError('signal_frequency must be a nonempty finite vector')
    idler = 2*frame-f
    if np.any(f <= 0) or np.any(idler <= 0):
        raise ValueError('signal and idler must have positive physical frequencies')
    matrix = _drift(r, detuning, pair_drive)
    eigenvalues, _, label = _stability(r, matrix)
    if label != 'stable':
        raise ValueError('small-signal response requires a strictly stable operating point')
    omega = 2*np.pi*(f-frame)
    inverse = np.linalg.inv(-1j*omega[:, None, None]*np.eye(2)-matrix)
    scattering = r.kappa_external*inverse-np.eye(2)
    bath = np.sqrt(r.kappa_external*r.kappa_internal)*inverse
    gain, conversion = abs(scattering[:, 0, 0])**2, abs(scattering[:, 0, 1])**2
    ns, ni = _occupation(f, r.internal_temperature), _occupation(idler, r.internal_temperature)
    numerator = .5*conversion+abs(bath[:, 0, 0])**2*(ns+.5)+abs(bath[:, 0, 1])**2*(ni+.5)
    noise = np.divide(numerator, gain, out=np.full_like(gain, np.inf), where=gain>0)
    degenerate = np.abs(f-frame) <= max(1e-6, abs(frame)*1e-14)
    noise[degenerate] = np.nan
    return JPAResponse(f, idler, frame, scattering, bath, gain, conversion,
                       noise, degenerate, eigenvalues, r)


def kerr_response(resonator, point, signal_frequency):
    """Signal/idler scattering about a selected Kerr pump operating point."""
    if point.resonator != resonator:
        raise ValueError('pump point belongs to a different resonator')
    delta = 2*np.pi*(point.frequency-resonator.frequency)-2*resonator.kerr*point.photons
    return _response(resonator, point.frequency, delta,
                     -1j*resonator.kerr*point.amplitude**2, signal_frequency)


def flux_response(resonator, pump_frequency, pair_drive, signal_frequency, *, frequency_shift=0.):
    """Effective flux-pumped DPO; frame is HALF the physical pump frequency.

    pair_drive is complex lambda [rad/s] in da/dt=...+lambda*a^dagger.
    frequency_shift is a supplied pump-induced resonance shift [Hz].
    Kerr saturation and flux-to-lambda circuit extraction are not included.
    """
    positive(pump_frequency, 'pump_frequency')
    if not np.isfinite(pair_drive) or not np.isfinite(frequency_shift):
        raise ValueError('pair_drive and frequency_shift must be finite')
    frame = pump_frequency/2
    delta = 2*np.pi*(frame-resonator.frequency-frequency_shift)
    return _response(resonator, frame, delta, pair_drive, signal_frequency)


@dataclass
class Quadratures:
    external_transfer: np.ndarray
    bath_transfer: np.ndarray
    output_covariance: np.ndarray  # Symmetrized x,p variance; vacuum = I/2.
    power_gains: np.ndarray  # Sorted singular values squared.
    noise_variances: np.ndarray  # Sorted covariance eigenvalues.


def quadratures(response):
    """DC quadrature response; supply exactly one degenerate frequency point.

    External input is vacuum; the internal bath uses resonator temperature.
    Output squeezing is spectral quadrature noise at zero analysis frequency,
    not intracavity variance. Away from degeneracy, use sideband scattering.
    """
    if len(response.signal_frequency) != 1 or not response.degenerate[0]:
        raise ValueError('quadratures requires a single degenerate frequency point')
    resonator = response.resonator
    transform = np.array([[1, 1j], [1, -1j]])/np.sqrt(2)
    ext = np.real_if_close(transform.conj().T@response.scattering[0]@transform).real
    bath = np.real_if_close(transform.conj().T@response.bath_scattering[0]@transform).real
    thermal = float(_occupation(np.array(response.frame_frequency), resonator.internal_temperature))
    covariance = .5*ext@ext.T+(thermal+.5)*bath@bath.T
    return Quadratures(ext, bath, covariance, np.sort(np.linalg.svd(ext, compute_uv=False)**2),
                       np.linalg.eigvalsh(covariance))


@dataclass
class JPATransient:
    time: np.ndarray
    amplitude: np.ndarray
    output: np.ndarray  # Complex output envelope, sqrt(photons/s).
    signal_gain: float
    idler_conversion: float


def kerr_transient(resonator, point, signal_frequency, signal_power, duration, *,
                   samples=4097, discard=.5, rtol=1e-8, atol=1e-10):
    """Classical pump+probe evolution and coherent sideband gain extraction.

    Starts at the selected stable pump state; discard selects settling time.
    Fits DC and +/- first three probe harmonics. Check duration/timestep and
    harmonic convergence for compression studies. No quantum noise is sampled.
    """
    if point.resonator != resonator or point.stability != 'stable':
        raise ValueError('a stable pump point for this resonator is required')
    positive(signal_frequency, 'signal_frequency'); positive(signal_power, 'signal_power')
    positive(duration, 'duration')
    if not isinstance(samples, int) or samples < 32 or not 0 <= discard < .9:
        raise ValueError('samples >= 32 and 0 <= discard < 0.9 required')
    r = resonator
    omega = 2*np.pi*(signal_frequency-point.frequency)
    if abs(omega)*duration*(1-discard) < 4*np.pi:
        raise ValueError('acquisition must contain at least two probe beat periods')
    if 2*point.frequency-signal_frequency <= 0:
        raise ValueError('idler frequency must be positive')
    if abs(omega)*duration/(samples-1)*3 >= np.pi:
        raise ValueError('sample interval aliases fitted probe harmonics')
    incoming = np.sqrt(point.power/(h*point.frequency))*np.exp(1j*point.phase)
    probe = np.sqrt(signal_power/(h*signal_frequency))
    delta = 2*np.pi*(point.frequency-r.frequency)
    scale = max(abs(point.amplitude), np.sqrt(r.kappa_external)*probe/r.kappa, 1.)
    def rhs(tau, state):
        a = state[0]*scale
        drive = incoming+probe*np.exp(-1j*omega*tau/r.kappa)
        rate = (1j*delta-r.kappa/2)*a-1j*r.kerr*abs(a)**2*a+np.sqrt(r.kappa_external)*drive
        return [rate/(scale*r.kappa)]
    time = np.linspace(0, duration, samples)
    sol = solve_ivp(rhs, (0., duration*r.kappa), [point.amplitude/scale],
                    t_eval=time*r.kappa, rtol=rtol, atol=atol,
                    max_step=min(.25, r.kappa/(12*abs(omega))))
    if not sol.success:
        raise RuntimeError(sol.message)
    amplitude = sol.y[0]*scale
    output = np.sqrt(r.kappa_external)*amplitude-incoming-probe*np.exp(-1j*omega*time)
    keep = time >= discard*duration
    harmonics = np.arange(-3, 4)
    fit = np.exp(-1j*omega*time[keep, None]*harmonics)
    coefficients = np.linalg.lstsq(fit, output[keep], rcond=None)[0]
    gain = abs(coefficients[4]/probe)**2
    conversion = abs(coefficients[2]/probe)**2
    return JPATransient(time, amplitude, output, float(gain), float(conversion))
