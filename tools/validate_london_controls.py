"""Native controls for kinetic energy, mutual sign and physical coil spacing."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import numpy as np

from twpa.em.fasthenry import FastHenry
from twpa.em.magnetic import jpa_magnetic_job
from twpa.layout import load_process
from twpa.layout.generate import JosephsonLayout,generate


def wire_job(scale):
    text=f'''* uniform thin-film London control
.units um
.default sigma=0
N0 x=0 y=0 z=0
N1 x=20 y=0 z=0
E0 N0 N1 w=2 h=.04 lambda={.45*scale:.16g} sigma=0 nwinc=1 nhinc=1 wx=0 wy=1 wz=0
.external N0 N1
.freq fmin=1e6 fmax=1e6 ndec=1
.end
'''
    return dict(schema_version=1,kind='uniform_london_wire',deck=text,
                source={'fixture':'20um x 2um x 40nm uniform-current wire'},manifest_sha256=None,
                stack={'lambda_um':.45*scale},settings={'frequencies_hz':[1e6]},
                ports=[dict(name='wire',positive='N0',negative='N1')],reduction=[[1]],patches=[],
                audit={'filaments':1},scope='Uniform-current analytic kinetic-energy control, not a device return-path model')


def run(root,installation=None,existing=False):
    root=Path(root)
    if not existing:root.mkdir(parents=True,exist_ok=False)
    solver=FastHenry.from_installation(installation);values={}
    for name,scale in [('base',1),('double',2),('pec',1e-4),('plus',np.sqrt(1.05)),('minus',np.sqrt(.95))]:
        value=json.loads((root/('wire-'+name)/'extraction.json').read_text()) if existing else solver.run(wire_job(scale),root/('wire-'+name))
        values[name]=value['samples'][0]['L_H'][0][0]
    expected=4e-7*np.pi*(.45e-6)**2*20e-6/(2e-6*.04e-6)
    derivative=(values['plus']-values['minus'])/.1
    ratio=(values['double']-values['pec'])/(values['base']-values['pec'])
    wire=dict(expected_kinetic_H=expected,derivative_kinetic_H=derivative,
              lambda_squared_ratio=ratio,passed=bool(abs(derivative/expected-1)<1e-8 and abs(ratio-4)<1e-6),values_H=values)
    records={}
    for name,gap in [('near',8.),('far',20.)]:
        gds=root/(name+'.gds')
        if not existing:generate(JosephsonLayout(inductive_coupler=True,coupler_gap=gap),load_process('jj_overlap_illustrative'),gds)
        job=jpa_magnetic_job(gds,pitch_um=4,local_pitch_um=1,nh=1,frequencies_hz=(1e6,))
        result=json.loads((root/name/'extraction.json').read_text()) if existing else solver.run(job,root/name)
        records[name]=np.array(result['samples'][0]['L_H'])
        if name=='near':
            reversed_job=deepcopy(job);port=reversed_job['ports'][2]
            old='.external '+port['positive']+' '+port['negative']
            port['positive'],port['negative']=port['negative'],port['positive']
            port['positive_patch'],port['negative_patch']=port['negative_patch'],port['positive_patch']
            reversed_job['deck']=reversed_job['deck'].replace(old,'.external '+port['positive']+' '+port['negative'])
            rev=json.loads((root/'reversed/extraction.json').read_text()) if existing else solver.run(reversed_job,root/'reversed')
            records['reversed']=np.array(rev['samples'][0]['L_H'])
    sign=np.diag([1.,-1.]);orientation=float(np.max(abs(records['reversed']-sign@records['near']@sign))/np.max(abs(records['near'])))
    coupling=dict(near_H=records['near'].tolist(),far_H=records['far'].tolist(),reversed_H=records['reversed'].tolist(),
                  orientation_relative_error=orientation,passed=bool(orientation<1e-10 and abs(records['far'][0,1])<abs(records['near'][0,1])))
    summary=dict(schema_version=1,wire=wire,coupling=coupling,simulator=solver.probe(),passed=wire['passed'] and coupling['passed'])
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2),flush=True)
    if not summary['passed']:raise RuntimeError('Native London/coupling controls failed')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--installation',type=Path)
    p.add_argument('--assess-existing',action='store_true',help='Assess preserved complete solver outputs without rerunning/overwriting them')
    a=p.parse_args();run(a.output,a.installation,a.assess_existing)
