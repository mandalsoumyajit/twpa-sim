"""Capture portable evidence from completed native London and JoSIM runs.

The reference is a compact report, not a sealed extraction to back-annotate
another GDS. Raw matrices and source hashes are preserved without host paths.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import numpy as np


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def compact(value):
    d=deepcopy(value)
    if 'artifact_sha256' in d:d['source_artifact_sha256']=d.pop('artifact_sha256')
    if 'run' in d:
        run=d['run'];run.pop('command',None);run['simulator'].pop('command',None)
    return d


def difference(a,b):
    result={}
    for key,field in [('signal_gain_db','gain'),('idler_conversion_db','idler_conversion')]:
        ratio=np.array(a['pumped']['measured'][field])/np.array(b['pumped']['measured'][field])
        result[key]=float(max(abs(10*np.log10(ratio))))
    return result


def capture(runs,output,figure):
    from twpa.em.coupled_jpa import CoupledJPA,pump,response
    m=read(runs/'magnetic-refined/magnetic-certified.json')
    e=read(runs/'magnetic-electric/electric-certified.json')
    cases={name:read(runs/f'{directory}/response.json') for name,directory in
           [('fine','magnetic-josim-fine'),('base','magnetic-josim-npy'),('window','magnetic-josim-window')]}
    final=cases['fine'];controls={}
    transient_runs={}
    for case,directory in [('fine','magnetic-josim-fine'),('base','magnetic-josim-npy'),('window','magnetic-josim-window')]:
        transient_runs[case]={}
        for name in ('passive-probe','pump-base','pump-probe'):
            root=runs/directory/name;manifest=read(root/'manifest.json');cmd=manifest['command']
            deck=(root/'input.cir').read_text(encoding='utf-8')
            if hashlib.sha256(deck.encode()).hexdigest()!=manifest['netlist_sha256']:
                raise ValueError('Native circuit input no longer matches its run manifest')
            transient_runs[case][name]=dict(analysis=int(cmd[cmd.index('-a')+1]),
                corrector=int(cmd[cmd.index('-c')+1]) if '-c' in cmd else 0,
                output_format=manifest['output_format'],netlist_sha256=manifest['netlist_sha256'],
                elapsed_seconds=manifest['elapsed_seconds'],returncode=manifest['returncode'],
                metadata=manifest['metadata'],netlist=deck)
    for name in m['certification']['comparison_artifacts']:
        controls[name]=compact(read(runs/f'magnetic-refined/magnetic-{name}/extraction.json'))
    coarse=read(runs/'magnetic-probe/primary/extraction.json')
    contact=read(runs/'magnetic-reference/magnetic-contact/extraction.json')
    a=np.array(coarse['samples'][0]['L_H']);b=np.array(contact['samples'][0]['L_H'])
    rejected=dict(primary_local_pitch_um=.5,comparison_local_pitch_um=1.,
                  max_self_relative=float(max(abs(np.diag(b-a)/np.diag(a)))),
                  mutual_relative=float(abs((b-a)[0,1]/a[0,1])),
                  self_gate=.02,mutual_gate=.05,passed=False,
                  source_artifacts=[coarse['artifact_sha256'],contact['artifact_sha256']])
    circuit=final['circuit'].copy();pec=np.array(m['pec_limit_H'])
    circuit.update(interconnect_inductance=pec[0,0],pickup_inductance=pec[1,1],mutual_inductance=pec[0,1])
    analytic=read(runs/'london-controls/summary.json');analytic['simulator'].pop('command',None)
    data=dict(schema_version=1,date='2026-10-28',magnetic=compact(m),magnetic_controls=controls,
              electric={k:e[k] for k in ('source','stack','capacitance_F','convergence')},
              electric_source_artifact_sha256=e['artifact_sha256'],response=final,
              refinement_cases=cases,transient_runs=transient_runs,sensitivity=dict(timestep=difference(cases['base'],final),
                  settling_window=difference(cases['base'],cases['window'])),
              pec_limit_modes=CoupledJPA(**circuit).modes(),analytic_and_coupling_controls=analytic,
              rejected_coarse_contact=rejected,
              galvanic_circuit_control=read(runs/'galvanic-control/response/response.json'),
              scope='Hybrid MIM electrostatic and full finite-conductor London MQS extraction; illustrative materials. '
                    'Galvanic case is a circuit-equation control with a coarse, uncertified magnetic mesh. '
                    'Not a full electric/package or fabricated-device validation.')
    if not all(v['passed'] for v in cases.values()) or not m['certification']['passed'] or not analytic['passed']:
        raise ValueError('Cannot capture failed reference')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if figure:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
        ax=axes[0,0];total=np.array(m['inductance_H'])*1e12;kin=np.array(m['kinetic_energy_H'])*1e12
        ax.bar(['Tank','Pickup'],np.diag(total-kin),label='Magnetic energy',color='#176b91')
        ax.bar(['Tank','Pickup'],np.diag(kin),bottom=np.diag(total-kin),label='Kinetic energy',color='#e6ae43')
        ax.set(ylabel='Inductance energy matrix (pH)',title=f'Signed mutual M = {total[0,1]:.3f} pH');ax.legend()
        ax=axes[0,1];ax.bar(['PEC limit','London'],[data['pec_limit_modes'][0]['frequency']/1e9,final['modes'][0]['frequency']/1e9],color=['#95a4ad','#176b91'])
        ax.set(ylabel='Lower mode (GHz)',title='Same circuit, no linewidth fitting')
        actual_circuit=CoupledJPA(**final['circuit'])
        ax=axes[1,0];p=final['passive'];f=np.array(p['measured']['signal_frequency'])/1e9
        dense=np.linspace(f.min(),f.max(),501);predicted=actual_circuit.passive(dense*1e9)
        for part,style in [('real','-'),('imag','--')]:
            ax.plot(dense,getattr(predicted,part),style,label='Model '+part)
            ax.plot(f,p['measured']['signal_'+part],'o',fillstyle='none',label='JoSIM '+part)
        ax.set(xlabel='Frequency (GHz)',ylabel='S11',title='Passive network');ax.legend(fontsize=8)
        ax=axes[1,1];p=final['pumped'];f=(np.array(p['measured']['signal_frequency'])-final['pump']['frequency'])/1e6
        dense=np.linspace(f.min(),f.max(),101)
        operating=pump(actual_circuit,final['pump']['frequency'],final['pump']['power'],harmonics=3)
        predicted=response(operating,final['pump']['frequency']+dense*1e6)
        for model,field,label in [('gain','gain','Signal'),('idler_conversion','idler_conversion','Idler')]:
            ax.plot(dense,10*np.log10(predicted[model]),label=label+' model')
            ax.plot(f,10*np.log10(p['measured'][field]),'o',fillstyle='none',label=label+' JoSIM')
        ax.set(xlabel='Probe offset (MHz)',ylabel='Gain / photon conversion (dB)',title='Full-sine pumped network');ax.legend(fontsize=8)
        for ax in axes.flat:ax.grid(alpha=.2)
        fig.suptitle('JPA: finite conductors → London matrix + openEMS → JoSIM',fontsize=13)
        figure.parent.mkdir(parents=True,exist_ok=True);fig.savefig(figure,dpi=150);plt.close(fig)
    print(json.dumps(dict(sensitivity=data['sensitivity'],rejected_contact=rejected,pec_modes=data['pec_limit_modes']),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,default=Path('validation/runs'))
    p.add_argument('--output',type=Path,default=Path('validation/magnetic-jpa-reference.json'))
    p.add_argument('--figure',type=Path,default=Path('docs/figures/magnetic-jpa-validation.png'))
    a=p.parse_args();capture(a.runs,a.output,a.figure)
