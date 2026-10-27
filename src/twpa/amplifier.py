"""Lossless analytic gain and distributed TLS/amplifier Langevin evolution."""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.constants import h, k as kb
from .linear import linear, supported, cascade
from .models import TLS


def default_frequency():
    return np.linspace(0, 11.5e9, 32769)


def select_pump(result, bias, spacing=1e9, offset=(2e6, 1.2e9), stopband_index=0):
    if not np.isfinite(spacing) or spacing < 0 or len(offset) != 2 or not 0 <= offset[0] < offset[1]:
        raise ValueError('invalid pump spacing or search offset')
    if not isinstance(stopband_index, int) or not 0 <= stopband_index < len(result.stopbands):
        raise ValueError('no selected stopband; supply pump_frequency explicitly')
    edge = result.stopbands[stopband_index, 1]
    f = result.frequency
    candidates = f[(f > edge+offset[0]) & (f < edge+offset[1]) & result.passband]
    fs, fi = (candidates-spacing)/2, (candidates+spacing)/2
    valid = supported(result, fs) & supported(result, fi)
    candidates, fs, fi = candidates[valid], fs[valid], fi[valid]
    if not len(candidates):
        raise ValueError('no supported pump candidates')
    k = PchipInterpolator(f, result.k, extrapolate=False)
    kp, ks, ki = k(candidates), k(fs), k(fi)
    mismatch = kp-ks-ki+bias.xi*bias.pump_current**2/8*(kp-2*ks-2*ki)
    return float(candidates[np.argmin(np.abs(mismatch))])


@dataclass
class GainResult:
    signal_frequency: np.ndarray
    pump_frequency: float
    gain: np.ndarray
    valid: np.ndarray

    @property
    def gain_db(self):
        return 10*np.log10(self.gain)


def small_signal(design, bias, signal_frequency, *, frequency=None, pump_frequency=None, **pump_options):
    r = linear(design, default_frequency() if frequency is None else frequency, bias.impedance)
    fp = select_pump(r, bias, **pump_options) if pump_frequency is None else pump_frequency
    fs = np.atleast_1d(np.asarray(signal_frequency, float))
    if fs.ndim != 1 or not len(fs) or np.any(~np.isfinite(fs)) or np.any(fs <= 0):
        raise ValueError('signal frequencies must be a nonempty positive finite vector')
    if not np.isfinite(fp) or not supported(r, fp):
        raise ValueError('pump must lie in a modeled passband')
    valid = supported(r, fs) & supported(r, fp-fs)
    interp = PchipInterpolator(r.frequency, r.k, extrapolate=False)
    kp, ks, ki = interp(fp), interp(fs[valid]), interp(fp-fs[valid])
    db = kp-ks-ki+bias.xi*bias.pump_current**2/8*(kp-2*ks-2*ki)
    g2 = (bias.epsilon*bias.pump_current/4)**2*ks*ki
    x = np.sqrt((g2-(db/2)**2).astype(complex))*design.total_cells
    sinhc = np.ones_like(x)
    nz = np.abs(x) > 1e-8
    sinhc[nz] = np.sinh(x[nz])/x[nz]
    sinhc[~nz] = 1+x[~nz]**2/6
    gain = np.full(fs.shape, np.nan)
    gain[valid] = np.abs(np.cosh(x)-1j*db/2*design.total_cells*sinhc)**2
    return GainResult(fs, fp, gain, valid)


@dataclass
class CoupledResult:
    signal_frequency: np.ndarray
    pump_frequency: float
    gain: np.ndarray
    added_noise: np.ndarray  # Input-referred symmetrized quanta; input idler vacuum.
    pump_output_current: np.ndarray
    transfer: np.ndarray
    bath_covariance: np.ndarray

    @property
    def gain_db(self):
        return 10*np.log10(self.gain)


def coupled(design, bias, signal_frequency, tls, *, frequency=None, pump_frequency=None,
            rtol=1e-7, atol=1e-10, pump_loss=True, **pump_options):
    """Local lossy Bloch parameters feed pump, signal/idler and bath-noise ODEs.

    Forward-wave, photon-normalized small-signal model; pump attenuates but is
    undepleted by mixing. Finite port reflections are not fed back into CMEs.
    The periodic medium is homogenized; tail cells contribute propagation length.
    """
    grid = default_frequency() if frequency is None else frequency
    ref = linear(design, grid, bias.impedance)
    fp = select_pump(ref, bias, **pump_options) if pump_frequency is None else pump_frequency
    fs = np.atleast_1d(np.asarray(signal_frequency, float))
    if fs.ndim != 1 or not len(fs) or np.any(~np.isfinite(fs)) or np.any(fs <= 0):
        raise ValueError('signal frequencies must be a nonempty positive finite vector')
    if not np.isfinite(fp) or not supported(ref, fp) or not np.all(supported(ref, fs) & supported(ref, fp-fs)):
        raise ValueError('pump, signal and idler must lie in modeled passbands')
    kinterp = PchipInterpolator(ref.frequency, ref.k)
    gains, noises, currents, transfers, covariances = [], [], [], [], []
    for signal in fs:
        modes = np.array([fp, signal, fp-signal])
        kref = kinterp(modes)

        def parameters(current):
            if tls.loss_tangent == 0 or not any(s.cell.participation for s in design.sections):
                return kref, np.zeros(3)
            t, scale = cascade(design.sections, modes, bias.impedance, tls=tls,
                               pump_current=current, pump_frequency=fp)
            trace = np.trace(t, axis1=1, axis2=2)*np.exp(scale)/2
            gamma = np.arccosh(trace.astype(complex))
            phase = gamma.imag+2*np.pi*np.round((kref*design.period_cells-gamma.imag)/(2*np.pi))
            kval, aval = phase/design.period_cells, gamma.real/design.period_cells
            if not pump_loss:
                kval[0], aval[0] = kref[0], 0.
            return kval, aval

        n = np.zeros(2) if tls.temperature == 0 else 1/np.expm1(np.minimum(h*modes[1:]/(kb*tls.temperature), 700))

        def rhs(position, state):
            current, phase = state[0].real, state[1].real
            k, alpha = parameters(max(0., current))
            kp, ks, ki = k
            mismatch = kp-ks-ki+bias.xi*current**2/8*(kp-2*ks-2*ki)
            coupling = bias.epsilon*current/4*np.sqrt(max(0., ks*ki))
            m = np.array([[-alpha[1], 1j*coupling*np.exp(1j*phase)],
                          [-1j*coupling*np.exp(-1j*phase), -alpha[2]]])
            u = state[2:6].reshape(2, 2)
            v = state[6:10].reshape(2, 2)
            diffusion = np.diag(2*alpha[1:]*(n+.5))
            return np.r_[-alpha[0]*current, mismatch, (m@u).ravel(),
                         (m@v+v@m.conj().T+diffusion).ravel()]

        initial = np.r_[bias.pump_current, 0., np.eye(2).ravel(), np.zeros(4)].astype(complex)
        sol = solve_ivp(rhs, (0., design.total_cells), initial, rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        last = sol.y[:, -1]
        u, v = last[2:6].reshape(2, 2), last[6:10].reshape(2, 2)
        gain = abs(u[0, 0])**2
        gains.append(gain)
        noises.append((.5*abs(u[0, 1])**2+v[0, 0].real)/gain)
        currents.append(last[0].real)
        transfers.append(u); covariances.append(v)
    return CoupledResult(fs, fp, np.array(gains), np.array(noises), np.array(currents),
                         np.array(transfers), np.array(covariances))
