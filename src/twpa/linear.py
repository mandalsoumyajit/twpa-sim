"""ABCD cascade, extended-zone Bloch dispersion, and finite-device S parameters."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import h, k as kb
from .models import TLS, positive


def tls_tangent(frequency, tls, pump_current=0., pump_frequency=0.):
    f = np.asarray(frequency, dtype=float)
    thermal = np.ones_like(f) if tls.temperature == 0 else np.tanh(h*f/(2*kb*tls.temperature))
    # A Lorentzian hole prevents a pump GHz away from saturating signal TLS.
    saturation = (tls.field_per_amp*pump_current/tls.critical_field)**2
    saturation = saturation/(1+((f-pump_frequency)/tls.linewidth)**2)
    return tls.loss_tangent*thermal/np.sqrt(1+saturation)


def cell_matrix(cell, frequency, impedance=50., tls=None, pump_current=0., pump_frequency=0.):
    positive(impedance, 'impedance')
    f = np.atleast_1d(np.asarray(frequency, dtype=float))
    if np.any(~np.isfinite(f)) or np.any(f < 0):
        raise ValueError('frequency must be finite and nonnegative')
    loss = 0 if tls is None else cell.participation*tls_tangent(f, tls, pump_current, pump_frequency)
    c = 2*cell.capacitance*(1-1j*loss)
    w = 2*np.pi*f
    den = 2-cell.stub_inductance*c*w*w
    if np.any(np.abs(den) < 1e-12):
        raise ValueError('frequency grid meets an effective-stub pole')
    t = np.zeros((len(f), 2, 2), complex)
    t[:, 0, 0] = 1
    t[:, 0, 1] = 1j*cell.inductance*w/impedance
    t[:, 1, 0] = 1j*2*c*w/den*impedance
    t[:, 1, 1] = 1-2*cell.inductance*c*w*w/den
    if cell.topology == 'symmetric':
        zy = t[:, 0, 1]*t[:, 1, 0]
        t[:, 0, 0] = 1+zy/2
        t[:, 1, 1] = 1+zy/2
        t[:, 0, 1] *= 1+zy/4
    return t


def _multiply(a, sa, b, sb):
    t = a @ b
    scale = np.max(np.abs(t), axis=(1, 2))
    if np.any(~np.isfinite(scale)) or np.any(scale == 0):
        raise ValueError('nonfinite or singular cascade')
    return t/scale[:, None, None], sa+sb+np.log(scale)


def _power(b, sb, count):
    a = np.broadcast_to(np.eye(2, dtype=complex), b.shape).copy()
    sa = np.zeros(len(b))
    while count:
        if count % 2:
            a, sa = _multiply(a, sa, b, sb)
        count //= 2
        if count:
            b, sb = _multiply(b, sb, b, sb)
    return a, sa


def cascade(sections, frequency, impedance=50., repeats=1, **kwargs):
    f = np.atleast_1d(frequency)
    t = np.broadcast_to(np.eye(2, dtype=complex), (len(f), 2, 2)).copy()
    scale = np.zeros(len(f))
    for s in sections:
        p, ps = _power(cell_matrix(s.cell, f, impedance, **kwargs), np.zeros(len(f)), s.count)
        t, scale = _multiply(t, scale, p, ps)
    return _power(t, scale, repeats)


@dataclass
class LinearResult:
    frequency: np.ndarray
    k: np.ndarray  # rad/cell, continuous extended-zone phase
    alpha: np.ndarray  # amplitude attenuation/cell
    passband: np.ndarray  # Classification from the corresponding lossless circuit
    stopbands: np.ndarray
    s11: np.ndarray
    s21: np.ndarray
    s21_db: np.ndarray
    total_cells: int


def linear(design, frequency, impedance=50., *, tls=None, pump_current=0., pump_frequency=0.):
    f = np.asarray(frequency, dtype=float)
    if f.ndim != 1 or len(f) < 2 or f[0] != 0 or np.any(~np.isfinite(f)) or np.any(np.diff(f) <= 0):
        raise ValueError('frequency must strictly increase from DC')
    base, bs = cascade(design.sections, f, impedance)
    base = base*np.exp(bs[:, None, None])
    ht = np.real(np.trace(base, axis1=1, axis2=2)/2)
    passband = np.abs(ht) <= 1+1e-12
    direction = np.sign(np.imag(base[:, 0, 1]+base[:, 1, 0]))
    direction[direction == 0] = 1
    q = np.unwrap(direction*np.arccos(np.clip(ht, -1, 1)))
    if np.any(np.diff(q) < -1e-7) or np.any(np.abs(np.diff(q)) > np.pi/2):
        raise ValueError('refine frequency grid or restrict range below stub resonances')
    kwargs = dict(tls=tls, pump_current=pump_current, pump_frequency=pump_frequency)
    if tls is not None and tls.loss_tangent and any(s.cell.participation for s in design.sections):
        t, scale = cascade(design.sections, f, impedance, **kwargs)
        trace = np.trace(t, axis1=1, axis2=2)*np.exp(scale)/2
        gamma = np.arccosh(trace.astype(complex))
        phase = np.unwrap(gamma.imag)
        # Anchor the extended-zone branch to the passive, lossless reference.
        phase += 2*np.pi*np.round((q-phase)/(2*np.pi))
        alpha = gamma.real/design.period_cells
        k = phase/design.period_cells
    else:
        k = q/design.period_cells
        alpha = np.arccosh(np.maximum(1, np.abs(ht)))/design.period_cells
    edges = np.diff(np.r_[False, ~passband, False].astype(int))
    bands = np.column_stack((f[np.where(edges == 1)[0]], f[np.where(edges == -1)[0]-1]))
    d, ds = cascade(design.sections, f, impedance, design.repeats, **kwargs)
    e, es = cascade(design.tail, f, impedance, **kwargs)
    d, ds = _multiply(d, ds, e, es)
    den = d.sum(axis=(1, 2))
    s21 = 2*np.exp(-ds)/den
    s11 = (d[:, 0, 0]+d[:, 0, 1]-d[:, 1, 0]-d[:, 1, 1])/den
    s21_db = 20/np.log(10)*(np.log(2)-ds-np.log(np.abs(den)))
    return LinearResult(f, k, alpha, passband, bands, s11, s21, s21_db, design.total_cells)


def supported(result, frequency):
    f = np.asarray(frequency)
    valid = (f > 0) & (f <= result.frequency[-1])
    for lo, hi in result.stopbands:
        valid &= ~((f >= lo) & (f <= hi))
    return valid


def network_scattering(blocks, frequency, impedance=50., **kwargs):
    """Cascade arbitrary (sections, repeat_count) blocks, e.g. rows and turns.

    Returns finite-network S11, S21, S21_dB without imposing periodicity.
    No single Bloch wavenumber or homogeneous amplifier gain is assigned.
    """
    f = np.atleast_1d(np.asarray(frequency, float))
    if f.ndim != 1 or np.any(~np.isfinite(f)) or np.any(f < 0):
        raise ValueError('frequency must be a finite nonnegative vector')
    d = np.broadcast_to(np.eye(2, dtype=complex), (len(f), 2, 2)).copy()
    ds = np.zeros(len(f))
    for sections, count in blocks:
        if not isinstance(count, int) or count < 0:
            raise ValueError('block repeat count must be a nonnegative integer')
        t, scale = cascade(sections, f, impedance, count, **kwargs)
        d, ds = _multiply(d, ds, t, scale)
    den = d.sum(axis=(1, 2))
    return dict(s11=(d[:, 0, 0]+d[:, 0, 1]-d[:, 1, 0]-d[:, 1, 1])/den,
                s21=2*np.exp(-ds)/den,
                s21_db=20/np.log(10)*(np.log(2)-ds-np.log(abs(den))))
