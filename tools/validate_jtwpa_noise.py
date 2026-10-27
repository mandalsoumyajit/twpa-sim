"""Record finite-circuit quantum-limit, loss-bath and spectral-resolution checks.

These are quantum consistency checks, not classical JoSIM noise simulations.
"""
import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import h, k as kb
from twpa.jpa_circuit import PHIBAR
from twpa.jtwpa import Chain, ResonantLoad, bloch, quantum_noise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args()
    base = Chain(PHIBAR/2.5e-9, 1e-12, 1e-12, 256)
    base = replace(base, port_impedance=float(bloch(base, [1e9]).impedance[0].real))
    rpm = replace(base, cells=384, resonator=ResonantLoad(1e-14, 1e-12,
                        1/((2*np.pi*1.008e9)**2*1e-12)))
    rpm = replace(rpm, port_impedance=float(bloch(rpm, [1e9]).impedance[0].real))
    loss = replace(base, junction_resistance=20000., ground_resistance=100000.)
    rpm_loss = replace(rpm, junction_resistance=20000., ground_resistance=100000.,
                        resonator=replace(rpm.resonator, resistance=100000.))
    cases = [('bare256', base, 0.), ('rpm384', rpm, 0.),
             ('rpm768', replace(rpm, cells=768), 0.),
             ('lossy256_cold', loss, 0.), ('lossy256_50mK', loss, .05),
             ('rpm_lossy384_50mK', rpm_loss, .05)]
    results = {}
    for name, chain, temperature in cases:
        power = .5*chain.port_impedance*(.4*chain.critical_current)**2
        options = dict(junction_temperature=temperature, ground_temperature=temperature,
                       resonator_temperature=temperature)
        r = quantum_noise(chain, 1e9, power, .9e9, **options)
        row = dict(chain=asdict(chain), pump_frequency=1e9, pump_power=power,
            signal_frequency=.9e9, temperatures=dict(port_temperature=0., **options),
            gain=r.gain, gain_db=r.gain_db, added_noise=r.added_noise,
            quantum_limit=r.quantum_limit, excess_noise=r.excess_noise,
            commutator_error=r.commutator_error, anomalous_commutator_error=r.anomalous_commutator_error,
            minimum_uncertainty_eigenvalue=r.minimum_uncertainty_eigenvalue,
            output_covariance_by_family={key:float(value[3, 3].real)
                                         for key, value in r.covariance_contributions.items()})
        assert r.phase_preserving and r.excess_noise > -1e-8
        assert max(r.commutator_error, r.anomalous_commutator_error) < 1e-8
        assert r.minimum_uncertainty_eigenvalue > -1e-8
        results[name] = row
        print(name, json.dumps(row), flush=True)
        if name in ('bare256', 'rpm384', 'rpm768'):
            refined = quantum_noise(chain, 1e9, power, .9e9, harmonics=5, quadrature=24)
            row['refined_added_noise'] = refined.added_noise
            row['added_noise_resolution_change'] = refined.added_noise-r.added_noise
            assert abs(row['added_noise_resolution_change']) < 1e-4
    assert results['lossy256_50mK']['added_noise'] > results['lossy256_cold']['added_noise']
    assert results['lossy256_cold']['excess_noise'] > results['bare256']['excess_noise']
    equilibrium = quantum_noise(rpm_loss, 1e9, 0., .9e9, port_temperature=.05,
            junction_temperature=.05, ground_temperature=.05, resonator_temperature=.05)
    occupation = np.repeat(1/np.expm1(h*equilibrium.frequency/(kb*.05))+.5, 2)
    equilibrium_error = float(np.max(abs(equilibrium.covariance-np.diag(occupation))))
    assert equilibrium_error < 1e-9
    payload = dict(schema_version=1, generated_utc=datetime.now(timezone.utc).isoformat(),
        kind='quantum consistency and convergence; not independent experimental data',
        model=dict(harmonics=3, sidebands=1, quadrature=16, nonlinearity='sine'),
        source_sha256=hashlib.sha256((Path(__file__).resolve().parents[1]/'src/twpa/jtwpa_noise.py').read_bytes()).hexdigest(),
        thermal_equilibrium_error=equilibrium_error, cases=results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2)+'\n')
    if args.plot:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        gains = np.geomspace(1, 1000, 200)
        axes[0].semilogx(gains, (1-1/gains)/2, label='Ideal phase-preserving limit')
        axes[0].axhline(.5, ls=':', color='gray', label='High-gain limit')
        for name in ['bare256', 'rpm384', 'rpm768']:
            row = results[name]
            axes[0].plot(row['gain'], row['added_noise'], 'o', label=name)
        axes[0].set_xlabel('Power gain G'); axes[0].set_ylabel('Input-referred added quanta')
        axes[0].legend(fontsize=8); axes[0].grid(alpha=.25)
        names = ['bare256', 'lossy256_cold', 'lossy256_50mK', 'rpm_lossy384_50mK']
        x = np.arange(len(names))
        axes[1].bar(x, [results[n]['added_noise'] for n in names], label='Added noise')
        axes[1].plot(x, [results[n]['quantum_limit'] for n in names], 'kx', ms=8, label='Limit at that gain')
        axes[1].set_xticks(x, ['Bare 256', 'Lossy 256\n0 K', 'Lossy 256\n50 mK', 'Lossy RPM\n50 mK'])
        axes[1].set_ylabel('Input-referred added quanta'); axes[1].legend(fontsize=8)
        axes[1].grid(axis='y', alpha=.25)
        fig.suptitle('4WM JTWPA quantum noise: finite gain and internal baths')
        args.plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.plot, dpi=180)


if __name__ == '__main__':
    main()
