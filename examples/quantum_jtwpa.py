"""Classical pump plus quantum vacuum/thermal noise in a lossy 4WM chain."""
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
