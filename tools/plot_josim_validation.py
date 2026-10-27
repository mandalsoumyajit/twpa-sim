"""Render the frozen circuit/model comparison; requires matplotlib."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from twpa import jpa
from twpa.jpa_circuit import ParallelJPA, SQUIDJPA

root = Path(__file__).resolve().parents[1]
data = json.loads((root/'tests/fixtures/josim_jpa.json').read_text())
rows = data['reference']['summary']['cases']
fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), layout='constrained')
blue, orange = '#2455a4', '#c45519'
for axis, key, title in zip(axes[0], ['kerr', 'flux'], ['Kerr-pumped JPA', 'Flux-pumped SQUID JPA']):
    row = rows[key]
    x = np.array([1, 2])
    dense = np.linspace(.75, 2.25, 151)
    if key == 'kerr':
        model = ParallelJPA(**row['circuit']).reduced_model()
        point = jpa.select_pump_point(jpa.pump_points(model, row['pump_frequency'], row['pump_power'], phase=np.pi/2))
        predicted = jpa.kerr_response(model, point, point.frequency+dense*1e6)
    else:
        circuit = SQUIDJPA(**row['circuit'])
        model = circuit.reduced_model()
        fp, amplitude, phase = row['flux_pump']
        predicted = jpa.flux_response(model, fp, circuit.pair_drive(amplitude, phase), fp/2+dense*1e6)
    axis.plot(dense, predicted.gain_db, '-', color=blue, label='Reduced model')
    axis.plot(x, 10*np.log10(row['measured']['gain']), 'x', color=orange, markersize=9, label='JoSIM circuit')
    axis.set(title=title, xlabel='Signal offset from frame (MHz)', ylabel='Signal gain (dB)', xticks=x)
    axis.legend(frameon=False)
power = [r['probe_power'] for r in rows['compression']['points']]
axes[1, 0].plot(10*np.log10(np.array(power)/1e-3),
                10*np.log10([r['model_gain'] for r in rows['compression']['points']]), '-o', color=blue)
axes[1, 0].plot(10*np.log10(np.array(power)/1e-3),
                10*np.log10([r['measured_gain'] for r in rows['compression']['points']]), 'x', color=orange, markersize=9)
axes[1, 0].set(title='Finite-signal Kerr response', xlabel='Incident probe power (dBm)', ylabel='Signal gain (dB)')
x = np.arange(2)
axes[1, 1].bar(x-.16, 10*np.log10(rows['flux']['model_quadrature_gains']), width=.32, color=blue, label='Reduced model')
axes[1, 1].bar(x+.16, 10*np.log10(rows['flux']['quadratures']['power_gains']), width=.32, color=orange, alpha=.8, label='JoSIM circuit')
axes[1, 1].set(title='Degenerate quadrature response', ylabel='Quadrature power gain (dB)',
               xticks=x, xticklabels=['Deamplified', 'Amplified'])
axes[1, 1].legend(frameon=False)
for axis in axes.ravel():
    axis.grid(alpha=.2)
fig.suptitle('Classical JPA validation: physical RCSJ circuits vs reduced models', fontsize=13)
out = root/'docs/figures/jpa-josim-validation.png'
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=170)
print(out)
