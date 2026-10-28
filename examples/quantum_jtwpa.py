"""Compare cold and warm internal loss baths in the same pumped junction line.

The finite 256-cell Josephson line is driven at 1 GHz and probed at 0.9 GHz,
with a 1.1 GHz four-wave-mixing idler. Junction and ground resistances cause
damping and must also supply fluctuations. Raising their temperature changes
thermal occupation, while circuit current/capacitance/resistance values stay
fixed. At zero temperature the baths still contribute vacuum fluctuations.

The output compares transmitted power gain in dB, input-referred added noise,
its finite-gain quantum lower bound, and excess above that bound. Noise is
reported in symmetrized photon units, not integrated output watts. This is
a Gaussian fluctuation calculation about a classical pump state.
Start with examples/howto/05_four_wave_jtwpa.py and docs/jtwpa-quantum-noise.md.
Physical limit: Caves et al., https://doi.org/10.1103/PhysRevA.86.063802.
"""
from twpa import jtwpa

chain = jtwpa.Chain(critical_current=1.3e-7, ground_capacitance=1e-12,
                    junction_capacitance=1e-12, cells=256,
                    junction_resistance=20000., ground_resistance=100000.)
for temperature in [0., .05]:
    result = jtwpa.quantum_noise(chain, 1e9, 7e-14, .9e9,
                    junction_temperature=temperature, ground_temperature=temperature)
    print(f'T={temperature:g} K, gain={result.gain_db:.4f} dB, '
          f'added={result.added_noise:.6f}, limit={result.quantum_limit:.6f}, '
          f'excess={result.excess_noise:.6f} quanta')
