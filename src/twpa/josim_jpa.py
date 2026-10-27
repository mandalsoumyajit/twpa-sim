"""Circuit measurement helpers for JoSIM resonator validation."""
from dataclasses import asdict
import numpy as np
from .jpa_circuit import Tone, jpa_netlist
from .measurement import fit_tones, port_waves


def run_circuit(runner, circuit, directory, tones, transient, *, flux_pump=None):
    metadata = dict(circuit_type=type(circuit).__name__, circuit=asdict(circuit),
                    tones=[asdict(t) for t in tones], transient=asdict(transient),
                    flux_pump=flux_pump)
    return runner.run(jpa_netlist(circuit, tones, transient, flux_pump=flux_pump),
                      directory, metadata=metadata)


def waveforms(result, impedance):
    voltage, drive = result.trace('V(PORT)'), result.trace('V(SRC0)')
    return port_waves(voltage, (drive-voltage)/impedance, impedance)


def compare_time_grids(result, baseline):
    if result.time.shape != baseline.time.shape or not np.allclose(result.time, baseline.time, rtol=0, atol=1e-18):
        raise ValueError('pump-only and probe runs must have identical time grids')


def sideband_measurement(result, baseline, impedance, frame_frequency, signal_frequencies):
    """Pump-only subtraction and simultaneous signal/idler coherent fitting.

    Returns power gains and photon-flux idler conversions. The incident waves
    are measured at the same reference plane, not assumed from source voltage.
    Additional mutually overlapping probes/idlers must be run separately.
    """
    compare_time_grids(result, baseline)
    fs = np.atleast_1d(np.asarray(signal_frequencies, float))
    fi = 2*frame_frequency-fs
    modes = np.r_[fs, fi]
    if np.any(fi <= 0) or len(np.unique(modes)) != len(modes):
        raise ValueError('probe and idler frequencies must be distinct and positive')
    frequencies = np.unique(np.r_[modes, frame_frequency, 2*frame_frequency, 3*frame_frequency])
    ain, aout = waveforms(result, impedance)
    bin_, bout = waveforms(baseline, impedance)
    incoming = fit_tones(result.time, ain-bin_, frequencies)
    outgoing = fit_tones(result.time, aout-bout, frequencies)
    si = np.searchsorted(frequencies, fs); ii = np.searchsorted(frequencies, fi)
    if np.any(abs(incoming[si]) <= 1e-20):
        raise ValueError('no measurable incident probe')
    signal = outgoing[si]/incoming[si]
    idler = outgoing[ii]/np.conj(incoming[si])*np.sqrt(fs/fi)
    return dict(signal_frequency=fs.tolist(), signal_real=signal.real.tolist(),
                signal_imag=signal.imag.tolist(), gain=(abs(signal)**2).tolist(),
                idler_real=idler.real.tolist(), idler_imag=idler.imag.tolist(),
                idler_conversion=(abs(idler)**2).tolist())


def degenerate_measurement(results, baseline, impedance, frequency):
    """Two orthogonal probe phases determine the real quadrature transfer."""
    if len(results) != 2:
        raise ValueError('two independently phased probe runs required')
    ins, outs = [], []
    bi, bo = waveforms(baseline, impedance)
    for result in results:
        compare_time_grids(result, baseline)
        ai, ao = waveforms(result, impedance)
        frequencies = [frequency, 2*frequency, 3*frequency]
        a = fit_tones(result.time, ai-bi, frequencies)[0]
        b = fit_tones(result.time, ao-bo, frequencies)[0]
        ins.append([a.real, a.imag]); outs.append([b.real, b.imag])
    incoming = np.array(ins).T
    if np.linalg.cond(incoming) > 1e6:
        raise ValueError('probe phases must span two independent quadratures')
    transfer = np.array(outs).T@np.linalg.inv(incoming)
    return dict(transfer=transfer.tolist(), power_gains=np.sort(np.linalg.svd(transfer, compute_uv=False)**2).tolist())
