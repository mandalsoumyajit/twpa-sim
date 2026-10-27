"""Check the physical SNAIL loop's small-geometric-inductance limit in JoSIM."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
from twpa import three_wave as tw
from twpa.backends import JoSIM
from twpa.jpa_circuit import Tone, Transient
from twpa.josim_jtwpa import measure_chain


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--executable',default='josim-cli')
    p.add_argument('--wsl')
    p.add_argument('--corrector',type=int,default=0)
    args=p.parse_args()
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    c=tw.example_chain('snail',16,dispersion_engineered=True)
    transient=Transient(5e-13,500e-9,400e-9,10e-12)
    runner=JoSIM(args.executable,args.wsl,timeout=600,corrector=args.corrector)
    frequencies=[.45e9,.55e9,1e9]
    exact=tw.linear_response(c,frequencies)
    rows=[]
    for name,inductance in [('nominal',1e-13),('half',5e-14)]:
        run=runner.run(tw.chain_netlist(c,[Tone(f,1e-22,ramp_time=30e-9) for f in frequencies],
                    transient,snail_loop_inductance=inductance),root/name,
                    metadata=dict(chain=asdict(c),snail_loop_inductance=inductance))
        m=measure_chain(run,c,frequencies)
        s21=m['transmitted']/m['incident']
        mean=float(np.mean(run.trace('P(BJS8)')))
        rows.append(dict(inductance=inductance,mean_phase=mean,
            equilibrium_error=mean-c.element.equilibrium,
            max_s21_error=float(np.max(abs(s21-exact['s21'])))))
    payload=dict(chain=asdict(c),simulator=runner.probe(),rows=rows)
    (root/'summary.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(payload,indent=2))
    if not abs(rows[1]['equilibrium_error'])<abs(rows[0]['equilibrium_error'])*.8:
        raise SystemExit('SNAIL small-loop equilibrium did not converge')


if __name__=='__main__':
    main()
