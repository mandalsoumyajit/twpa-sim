"""Two-port coherent measurements for physical JoSIM JTWPA chains."""
from dataclasses import asdict
import numpy as np
from .jtwpa import chain_netlist
from .measurement import fit_tones, port_waves


def run_chain(runner, chain, directory, tones, transient):
    return runner.run(chain_netlist(chain, tones, transient), directory,
                      metadata=dict(chain=asdict(chain), tones=[asdict(t) for t in tones],
                                    transient=asdict(transient)))


def measure_chain(result, chain, frequencies, baseline=None):
    """Return physical two-port phasors; output load has no incident wave."""
    f = np.atleast_1d(np.asarray(frequencies, float))
    time = result.time
    drive = result.trace('V(SRC0)')
    vin = result.trace('V(N0)')
    vout = result.trace(f'V(N{chain.cells})')
    if baseline is not None:
        if baseline.time.shape != time.shape or not np.allclose(baseline.time, time, atol=1e-18, rtol=0):
            raise ValueError('baseline and probe time grids differ')
        drive = drive-baseline.trace('V(SRC0)')
        vin = vin-baseline.trace('V(N0)')
        vout = vout-baseline.trace(f'V(N{chain.cells})')
    incident, reflected = port_waves(vin, (drive-vin)/chain.port_impedance, chain.port_impedance)
    return dict(frequency=f, incident=fit_tones(time, incident, f),
                reflected=fit_tones(time, reflected, f), transmitted=fit_tones(time, vout, f))


def measure_bloch(result, chain, frequencies):
    """Extract trace from three adjacent interior node voltages, including reflections."""
    if chain.cells < 4:
        raise ValueError('Bloch extraction needs interior neighbors')
    m = max(1, chain.cells//2)
    fields = [fit_tones(result.time, result.trace(f'V(N{j})'), frequencies) for j in (m-1, m, m+1)]
    trace = (fields[0]+fields[2])/(2*fields[1])
    return np.arccos(trace.astype(complex))
