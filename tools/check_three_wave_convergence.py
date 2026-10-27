"""Refit actual acquisition half-windows and check depleted spectral resolution."""
import argparse
import json
from pathlib import Path
import numpy as np
from twpa import three_wave as tw
from twpa.jtwpa import ResonantLoad
from twpa.backends.josim import load_csv, SimulationResult
from twpa.josim_jtwpa import measure_chain


def load_chain(summary):
    data=dict(summary['chain'])
    cls=tw.RFSQUID if summary['kind']=='rf-squid' else tw.SNAIL
    data['element']=cls(**data['element'])
    if data['resonator'] is not None:
        data['resonator']=ResonantLoad(**data['resonator'])
    return tw.ThreeWaveChain(**data)


def windows(root):
    d=json.loads((root/'summary.json').read_text());c=load_chain(d)
    actual=load_csv(root/'probe/output.csv');base=load_csv(root/d.get('baseline_directory','pump')/'output.csv')
    time=next(iter(actual.values()));middle=(time[0]+time[-1])/2
    fp,fs=d['pump_frequency'],d['signal_frequency'];fi=fp-fs
    result=[]
    for label,selection in [('early',time<=middle),('late',time>middle)]:
        def record(data):
            return SimulationResult(root,{k:v[selection] for k,v in data.items()},{})
        a,b=record(actual),record(base)
        m=measure_chain(a,c,[fs,fi,fp,2*fp,3*fp],b)
        pump=measure_chain(a,c,[fp,2*fp,3*fp])
        result.append(dict(window=label,gain=float(abs(m['transmitted'][0]/m['incident'][0])**2),
            conversion=float(abs(m['transmitted'][1]/m['incident'][0])**2*fs/fi),
            pump=float(abs(pump['transmitted'][0]/pump['incident'][0])**2)))
    return dict(windows=result,late_minus_early_db={key:float(10*np.log10(result[1][key]/result[0][key]))
                                                  for key in ['gain','conversion','pump']})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    out=dict(windows={name:windows(args.root/name) for name in
                      ['rf128-fine','snail128-fine','rf128-loss','snail128-loss-depleted']})
    summary=json.loads((args.root/'snail128-loss-depleted/summary.json').read_text())
    c=load_chain(summary);a=summary['active']
    rows=[]
    for h,b,q in [(4,3,32),(3,4,32)]:
        r=tw.three_wave_response(c,summary['pump_frequency'],summary['pump_power'],
            summary['signal_frequency'],summary['signal_power'],depleted=True,
            harmonics=h,sidebands=b,quadrature=q)
        rows.append(dict(harmonics=h,sidebands=b,quadrature=q,
            gain_change_db=float(10*np.log10(r.gain/a['model_gain'])),
            conversion_change_db=float(10*np.log10(r.idler_conversion/a['model_conversion'])),
            pump_change_db=float(10*np.log10(r.pump_transmission/a['model_pump_transmission']))))
    out['depleted_resolution']=rows
    out['timestep_settling_refinement']={}
    for kind in ['rf','snail']:
        fine=json.loads((args.root/f'{kind}128-fine/summary.json').read_text())['active']
        coarse=json.loads((args.root/f'{kind}128-smooth/summary.json').read_text())['active']
        out['timestep_settling_refinement'][kind]={key:float(10*np.log10(fine[field]/coarse[field])) for key,field in
            [('gain_db','measured_gain'),('conversion_db','measured_conversion'),('pump_db','measured_pump_transmission')]}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
