"""3WM rf-SQUID/SNAIL quantum-limit, loss and convergence benchmarks."""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import h, k as kb
from twpa import three_wave as tw


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--plot', type=Path)
    args = p.parse_args()
    cases = {}
    for kind in ['rf-squid', 'snail']:
        for name, cells, resistance, ground, temperature in [
            ('128-cold',128,1e9,None,0.), ('384-cold',384,1e9,None,0.),
            ('128-loss-cold',128,20000.,100000.,0.), ('128-loss-hot',128,20000.,100000.,.05)]:
            c = tw.example_chain(kind, cells, junction_resistance=resistance,
                                 ground_resistance=ground, dispersion_engineered=True)
            power = .5*c.port_impedance*(.18*c.linear_chain.critical_current)**2
            q = tw.quantum_noise(c, 1e9, power, .45e9,
                                junction_temperature=temperature, ground_temperature=temperature)
            row = dict(kind=kind, chain=asdict(c), pump_power=power, pump_frequency=1e9,
                signal_frequency=.45e9, internal_temperature=temperature, gain=q.gain,
                gain_db=q.gain_db, added_noise=q.added_noise, quantum_limit=q.quantum_limit,
                excess_noise=q.excess_noise, commutator_error=q.commutator_error,
                anomalous_commutator_error=q.anomalous_commutator_error,
                minimum_uncertainty_eigenvalue=q.minimum_uncertainty_eigenvalue,
                output_covariance_by_family={key:float(value[1,1].real) for key,value in q.covariance_contributions.items()})
            assert q.phase_preserving and q.excess_noise > -1e-8
            assert max(q.commutator_error,q.anomalous_commutator_error)<1e-8
            assert q.minimum_uncertainty_eigenvalue>-1e-8
            if name=='384-cold':
                assert q.gain>100 and 0<=q.excess_noise<.02 and abs(q.added_noise-.5)<.02
            if name=='128-cold':
                refined = tw.quantum_noise(c, 1e9, power, .45e9, harmonics=4, quadrature=32)
                row['gain_resolution_change_db'] = float(10*np.log10(refined.gain/q.gain))
                row['noise_resolution_change'] = refined.added_noise-q.added_noise
                assert abs(row['gain_resolution_change_db'])<.03
                assert abs(row['noise_resolution_change'])<.001
                # Common-temperature equilibrium is an independent passive FDT check.
                eq = tw.quantum_noise(c,1e9,0.,.45e9,port_temperature=.05,
                                     junction_temperature=.05,resonator_temperature=.05)
                expected = np.repeat(1/np.expm1(h*eq.frequency/(kb*.05))+.5,2)
                row['equilibrium_error'] = float(np.max(abs(eq.covariance-np.diag(expected))))
                assert row['equilibrium_error']<1e-9
            cases[f'{kind}-{name}'] = row
            print(f'{kind}-{name}',json.dumps(row),flush=True)
        assert cases[f'{kind}-128-loss-hot']['added_noise']>cases[f'{kind}-128-loss-cold']['added_noise']
        assert cases[f'{kind}-128-loss-cold']['excess_noise']>cases[f'{kind}-128-cold']['excess_noise']
    root=Path(__file__).resolve().parents[1]
    payload=dict(kind='quantum consistency calculations, not independent measurements',
        generated_utc=datetime.now(timezone.utc).isoformat(),
        controls=dict(harmonics=3,sidebands=1,quadrature=24),
        source_sha256={file:hashlib.sha256((root/file).read_bytes()).hexdigest()
            for file in ['src/twpa/three_wave.py','src/twpa/jtwpa_noise.py']},cases=cases)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2)+'\n')
    if args.plot:
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
        gains=np.geomspace(1,1000,300)
        axes[0].semilogx(gains,(1-1/gains)/2,'k-',label='Finite-gain quantum limit')
        for kind,marker in [('rf-squid','o'),('snail','x')]:
            rows=[cases[f'{kind}-{name}'] for name in ['128-cold','384-cold']]
            axes[0].semilogx([r['gain'] for r in rows],[r['added_noise'] for r in rows],marker,label=kind)
        axes[0].set_xlabel('Signal power gain G');axes[0].set_ylabel('Input-referred added quanta')
        axes[0].legend();axes[0].grid(alpha=.25)
        names=['128-cold','128-loss-cold','128-loss-hot'];x=np.arange(3)
        for i,kind in enumerate(['rf-squid','snail']):
            axes[1].bar(x+(i-.5)*.35,[cases[f'{kind}-{name}']['added_noise'] for name in names],.35,label=kind)
        axes[1].set_xticks(x,['Small loss, 0 K','Internal loss, 0 K','Internal loss, 50 mK'])
        axes[1].set_ylabel('Added quanta, 128 cells');axes[1].legend();axes[1].grid(axis='y',alpha=.25)
        args.plot.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.plot,dpi=180)


if __name__=='__main__':
    main()
