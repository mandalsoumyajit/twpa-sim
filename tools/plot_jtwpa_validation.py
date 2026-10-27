"""Render measured JoSIM/model comparisons from the frozen fixture."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixture', type=Path, default=Path('tests/fixtures/josim_jtwpa.json'))
    p.add_argument('--output', type=Path, default=Path('docs/jtwpa-validation.png'))
    args = p.parse_args()
    data = json.loads(args.fixture.read_text())['cases']
    names = ['bare256', 'bare384', 'depleted256', 'rpm384']
    labels = ['Bare 256', 'Bare 384', 'Depleted 256', 'RPM 384']
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    x = np.arange(len(names))
    for ax, field, ylabel in [(axes[0], 'gain', 'Signal power gain (dB)'),
                              (axes[1], 'conversion', 'Idler photon conversion (dB)')]:
        for prefix, marker in [('model', 'o'), ('measured', 'x')]:
            values = [data[n]['summary']['active'][f'{prefix}_{field}'] for n in names]
            ax.plot(x, 10*np.log10(values), marker, ms=8, label='Model' if prefix=='model' else 'JoSIM')
        ax.set_xticks(x, labels, rotation=25, ha='right')
        ax.set_ylabel(ylabel); ax.grid(alpha=.25); ax.legend()
    d = data['depleted256']['summary']['active']
    f = np.array(d['retained_frequency'])/1e9
    order = np.argsort(f)
    for prefix, marker in [('model', 'o'), ('measured', 'x')]:
        axes[2].plot(f[order], 10*np.log10(np.array(d[f'{prefix}_output_power'])[order]/1e-3),
                     marker, ms=7, label='Model' if prefix=='model' else 'JoSIM')
    axes[2].set_xlabel('Frequency (GHz)'); axes[2].set_ylabel('Depleted output tone power (dBm)')
    axes[2].grid(alpha=.25); axes[2].legend()
    fig.suptitle('4WM Josephson chains: finite-circuit model vs RCSJ transients')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)


if __name__ == '__main__':
    main()
