"""Pump-depleted three-wave CMEs ported from sim_alvin_twpa.py.

Spatial phase is integrated as in sim_gradients_exact.py, avoiding hardcoded
routing and using user-supplied normalized-coordinate profiles.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from .amplifier import default_frequency, select_pump
from .linear import linear, supported
from .models import positive


def cme_rhs(position, currents, epsilon, wavenumbers, phase=None):
    """Original nonlinear RHS, with optional accumulated spatial WKB phase."""
    kp, ks, ki = wavenumbers
    e1, e2 = epsilon
    ip, isig, iid = currents
    p2, s2, i2 = abs(ip)**2, abs(isig)**2, abs(iid)**2
    phi = (kp-ks-ki)*position if phase is None else phase
    return np.array([
        1j*kp*e1*isig*iid*np.exp(-1j*phi)/4+1j*kp*e2*ip*(p2+2*s2+2*i2)/8,
        1j*ks*e1*ip*np.conj(iid)*np.exp(1j*phi)/4+1j*ks*e2*isig*(2*p2+s2+2*i2)/8,
        1j*ki*e1*ip*np.conj(isig)*np.exp(1j*phi)/4+1j*ki*e2*iid*(2*p2+2*s2+i2)/8])


@dataclass
class NonlinearResult:
    position: np.ndarray
    currents: np.ndarray  # (pump, signal, idler) x position; peak amperes
    gain: float
    pump_frequency: float


def route_profile(mode='ideal', cycles=1):
    """Return u(s), s in [0,1]. Modes map to original exact-gradient cases."""
    positive(cycles, 'cycles')
    profiles = {'ideal': lambda s: 0., 'linear': lambda s: s-.5,
                'from-input': lambda s: s,
                'sine': lambda s: .5*np.sin(2*np.pi*cycles*s),
                'folded': lambda s: .5*np.cos(2*np.pi*cycles*s),
                'spiral': lambda s: .5*np.sqrt(s)*np.sin(2*np.pi*cycles*np.sqrt(s))}
    if mode not in profiles:
        raise ValueError(f'Unknown route mode: {mode}')
    return profiles[mode]


def nonlinear(design, bias, signal_frequency, signal_current, *, idler_current=0j,
              frequency=None, pump_frequency=None, gradient=0., profile=None,
              derivative_step=.001, samples=201, rtol=1e-8, atol=1e-12):
    """Integrate full pump/signal/idler CMEs, optionally with a film gradient.

    gradient is fractional inductance variation times profile(cell/Nc).
    dk/d(fractional L) is a centered finite difference scaling all circuit
    inductances. Nonlinear coefficients stay fixed, as in the source scripts.
    Currents are internally normalized to the largest input amplitude so the
    absolute tolerance does not overwhelm a small probe current.
    """
    positive(abs(signal_current), 'signal_current')
    if not np.isfinite(gradient) or not np.isfinite(idler_current):
        raise ValueError('gradient and idler current must be finite')
    positive(derivative_step, 'derivative_step')
    if derivative_step >= 1 or not isinstance(samples, int) or samples < 2:
        raise ValueError('derivative_step must be < 1 and samples >= 2')
    grid = default_frequency() if frequency is None else frequency
    ref = linear(design, grid, bias.impedance)
    fp = select_pump(ref, bias) if pump_frequency is None else pump_frequency
    modes = np.array([fp, signal_frequency, fp-signal_frequency])
    if not np.all(np.isfinite(modes)) or not np.all(supported(ref, modes)):
        raise ValueError('pump, signal and idler must lie in modeled passbands')
    kv = PchipInterpolator(ref.frequency, ref.k)(modes)
    derivative = np.zeros(3)
    if gradient:
        def shifted(factor):
            def shift(s):
                return replace(s, cell=replace(s.cell, inductance=s.cell.inductance*factor,
                               stub_inductance=s.cell.stub_inductance*factor))
            return replace(design, sections=tuple(shift(s) for s in design.sections),
                           tail=tuple(shift(s) for s in design.tail))
        plus = linear(shifted(1+derivative_step), grid, bias.impedance)
        minus = linear(shifted(1-derivative_step), grid, bias.impedance)
        if not np.all(supported(plus, modes) & supported(minus, modes)):
            raise ValueError('gradient derivative crosses a stopband')
        derivative = (PchipInterpolator(grid, plus.k)(modes)-PchipInterpolator(grid, minus.k)(modes))/(2*derivative_step)
    profile = route_profile('linear') if profile is None else profile
    current_scale = max(bias.pump_current, abs(signal_current), abs(idler_current))
    def rhs(t, y):
        kt = kv+derivative*gradient*profile(t/design.total_cells)
        rates = cme_rhs(t, y[:3]*current_scale, (bias.epsilon, bias.xi), kt, y[3].real)
        return np.r_[rates/current_scale, kt[0]-kt[1]-kt[2]]
    initial = np.array([bias.pump_current/current_scale, signal_current/current_scale,
                        idler_current/current_scale, 0.], complex)
    pos = np.linspace(0, design.total_cells, samples)
    sol = solve_ivp(rhs, (0, design.total_cells), initial, t_eval=pos, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    currents = sol.y[:3]*current_scale
    return NonlinearResult(sol.t, currents, float(abs(currents[1, -1]/signal_current)**2), fp)


def compression(design, bias, signal_frequency, input_power_dbm, **options):
    """Power sweep and first interpolated 1-dB compression crossing; NaN if absent.

    Small-signal reference is a separate low-probe-current integration.
    """
    powers = np.asarray(input_power_dbm, float)
    if powers.ndim != 1 or len(powers)<2 or np.any(~np.isfinite(powers)) or np.any(np.diff(powers)<=0):
        raise ValueError('input powers must strictly increase')
    reference = nonlinear(design, bias, signal_frequency, max(bias.pump_current*1e-6, 1e-14), **options)
    gains = np.array([nonlinear(design, bias, signal_frequency,
                      np.sqrt(2*1e-3*10**(p/10)/bias.impedance), **options).gain for p in powers])
    deficit = 10*np.log10(reference.gain/gains)
    crossing = np.flatnonzero(deficit >= 1)
    p1db = np.nan
    if len(crossing) and crossing[0] > 0:
        j = crossing[0]
        p1db = powers[j-1]+(powers[j]-powers[j-1])*(1-deficit[j-1])/(deficit[j]-deficit[j-1])
    return gains, float(p1db)
