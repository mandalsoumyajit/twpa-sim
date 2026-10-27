"""Coherent real-waveform tone fitting and power-wave reference planes."""
import numpy as np


def fit_tones(time, values, frequencies):
    """Peak phasors with exp(-i omega t): real(v)=Re[V exp(-i omega t)].

    DC and all supplied tones are fitted simultaneously. An irregular time
    grid is allowed. Frequencies must be resolved by the acquisition duration.
    """
    t, y, f = np.asarray(time, float), np.asarray(values, float), np.asarray(frequencies, float)
    if t.ndim != 1 or y.shape != t.shape or f.ndim != 1 or not len(f):
        raise ValueError('time, waveform and frequency vectors required')
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(y)) or not np.all(np.isfinite(f)):
        raise ValueError('tone fit inputs must be finite')
    if len(t) < 2*len(f)+2 or np.any(np.diff(t) <= 0) or np.any(f <= 0) or len(np.unique(f)) != len(f):
        raise ValueError('insufficient samples, invalid time or duplicate/nonpositive tones')
    if np.max(np.diff(t))*max(f) >= .5:
        raise ValueError('sampling interval aliases fitted tones')
    span = t[-1]-t[0]
    if len(f)>1 and min(np.diff(np.sort(f)))*span < 1.:
        raise ValueError('acquisition does not resolve separate tones')
    angles = 2*np.pi*t[:, None]*f
    design = np.column_stack([np.ones(len(t)), np.cos(angles), np.sin(angles)])
    solution, _, rank, singular = np.linalg.lstsq(design, y, rcond=None)
    if rank < design.shape[1] or singular[0]/singular[-1] > 1e8:
        raise ValueError('ill-conditioned tone fit')
    phasors = solution[1:1+len(f)]+1j*solution[1+len(f):]
    return phasors


def port_waves(voltage, current, impedance):
    """Peak voltage waves, current positive INTO the device.

    Incident=(V+Z0 I)/2, reflected=(V-Z0 I)/2; power is |wave|^2/(2 Z0).
    The same expression applies to instantaneous traces or fitted phasors.
    """
    if not np.isfinite(impedance) or impedance <= 0:
        raise ValueError('positive finite reference impedance required')
    return ((np.asarray(voltage)+impedance*np.asarray(current))/2,
            (np.asarray(voltage)-impedance*np.asarray(current))/2)
