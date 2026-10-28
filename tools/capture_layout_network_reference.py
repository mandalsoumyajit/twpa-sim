"""Capture portable evidence from the completed multi-family benchmark tree.

The report preserves native matrices, input decks and source hashes without
host paths. It is not a transferable extraction certificate for another mask.
"""
import argparse
from datetime import date
import hashlib
import json
from pathlib import Path
import numpy as np
from twpa.em import LayoutNetwork,load_circuit,certify_magnetic
from twpa.em.magnetic_annotation import read_magnetic
from twpa.em.backannotation import read_extraction
from twpa.layout.generate import digest


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def portable(value):
    if isinstance(value,list):return [portable(v) for v in value]
    if not isinstance(value,dict):return value
    return {('source_artifact_sha256' if k=='artifact_sha256' else k):portable(v)
            for k,v in value.items() if k not in ('command','report')}


def magnetic_difference(a,b):
    a=np.array(a['samples'][0]['L_H']);b=np.array(b['samples'][0]['L_H'])
    scale=np.sqrt(np.outer(np.diag(a),np.diag(a)));off=~np.eye(len(a),dtype=bool)
    strong=off&(abs(a)/scale>=.01);delta=b-a
    return dict(max_self_relative=float(max(abs(np.diag(delta)/np.diag(a)))),
        mutual_relative=float(max(abs(delta[strong]/a[strong]))) if np.any(strong) else 0.,
        coupling_absolute=float(max(abs(delta[off])/scale[off])))


def case(runs,name,field_directory,primary_directory):
    folder=runs/(name+'-accepted');gds=runs/(name+'.gds')
    manifest=read(gds.with_suffix('.manifest.json'))
    ann=read(folder/'backannotation.json');circuit=load_circuit(folder/'backannotation.json')
    if not isinstance(circuit,LayoutNetwork):raise ValueError('Expected terminal network')
    if ann['source']['gds_sha256']!=sha(gds):raise ValueError('Layout source changed')
    magnetic=read_magnetic(runs/field_directory/'certified.json')
    controls={k:read_magnetic(runs/field_directory/k/'extraction.json') for k in ('bulk','contact','thickness','frequency')}
    primary=read_magnetic(runs/primary_directory/'extraction.json')
    if certify_magnetic(primary,controls)['artifact_sha256']!=magnetic['artifact_sha256'] or not magnetic['certification']['passed']:
        raise ValueError('Native field evidence does not reproduce the certificate')
    if ann['magnetic_sha256']!=magnetic['artifact_sha256']:raise ValueError('Different annotated magnetic matrix')
    electric={c['id']:read_extraction(folder/('electric-'+c['id']+'.json')) for c in manifest['components'] if c['kind']=='CAP'}
    if {k:v['artifact_sha256'] for k,v in electric.items()}!=ann['electric_sha256']:raise ValueError('Different capacitor evidence')
    responses={};native={}
    for label in ('primary','refined','flux-pump'):
        path=folder/label/'response.json'
        if not path.exists():
            if label=='flux-pump':continue
            raise FileNotFoundError(path)
        d=read(path)
        if not d['passed'] or d['circuit']!=ann['circuit']:raise ValueError('Failed/different circuit response')
        native[label]={}
        for simulation in ('passive','pump-base','pump-probe'):
            root=folder/label/simulation;m=read(root/'manifest.json');deck=(root/'input.cir').read_text(encoding='utf-8')
            if sha(root/'input.cir')!=m['netlist_sha256'] or m['returncode']!=0 or m['metadata']['circuit']!=ann['circuit']:
                raise ValueError('Native netlist/provenance differs')
            command=m['command']
            native[label][simulation]=dict(netlist=deck,netlist_sha256=m['netlist_sha256'],
                analysis=int(command[command.index('-a')+1]),corrector=int(command[command.index('-c')+1]) if '-c' in command else 0,
                output_format=m['output_format'],returncode=m['returncode'],elapsed_seconds=m['elapsed_seconds'],
                started_utc=m['started_utc'],source_manifest_sha256=sha(root/'manifest.json'))
        d.pop('circuit');d.pop('native_runs',None);d['source_response_sha256']=sha(path);responses[label]=portable(d)
    step=read(folder/'timestep.json')
    for k,error in step['changes_db'].items():
        actual=float(10*np.log10(responses['refined']['pumped'][k]/responses['primary']['pumped'][k]))
        if not np.isclose(error,actual,atol=1e-12):raise ValueError('Stale timestep comparison')
    if not step['passed']:raise ValueError('Timestep check failed')
    drc=read(runs/'native-drc'/(name+'.json'))
    if not drc['passed'] or drc['violations']:raise ValueError('Native DRC failed')
    return dict(design=manifest['design'],family=manifest['family'],process_id=manifest['process']['id'],
        process_status=manifest['process']['status'],source=ann['source'],manifest_sha256=sha(gds.with_suffix('.manifest.json')),
        magnetic=portable(magnetic),magnetic_controls={k:portable(v) for k,v in controls.items()},
        electric={k:portable(v) for k,v in electric.items()},backannotation=portable(ann),responses=responses,
        timestep=step,transient_runs=native,native_drc=dict(passed=True,violations=0,
        report_sha256=sha(runs/'native-drc'/(name+'.lyrdb')),gds_sha256=sha(gds)))


def capture(runs,kit,output):
    configurations=[('squid_jpa','squid-patch1-fine-controls','squid-patch1-fine'),
        ('jtwpa4wm','four-patch1-fine-controls','four-patch1-fine'),
        ('rfsquid3wm','rfsquid3wm-final-fields','rfsquid3wm-final-fields/primary'),
        ('snail-two-narrow','snail-two-narrow-final-fields','snail-two-narrow-final-fields/primary')]
    cases={name:case(runs,name,fields,primary) for name,fields,primary in configurations}
    # Preserve rejected native attempts rather than changing their gates.
    rejected=[portable(read_magnetic(runs/p)) for p in
        ('squid-controls-v3/magnetic-certified.json','squid-patch1-controls/certified.json')]
    if any(d['certification']['passed'] for d in rejected):raise ValueError('Expected rejected coarse controls')
    a=read_magnetic(runs/'snail-two-compact-final-fields/primary/extraction.json')
    b=read_magnetic(runs/'snail-two-compact-final-fields/contact/extraction.json')
    if a['source']!=b['source'] or a['patches']!=b['patches']:raise ValueError('Rejected contact comparison changes the physical problem')
    snail_rejected=dict(reason='Fixed-terminal contact sensitivity exceeds the unchanged normalized-mutual gate',
        design=read(runs/'snail-two-compact.manifest.json')['design'],
        primary=portable(a),comparison=portable(b),checks=magnetic_difference(a,b),coupling_gate=.005,passed=False)
    failed_electric={}
    for name,path in [('one_ohm_low_frequency','snail-two-compact-broadband-electric/electric-CAP4-primary/extraction.json'),
                      ('short_ten_ohm_transient','snail-r10-broadband-controls/duration/extraction.json')]:
        p=runs/path;d=read(p);d['run']=read(p.parent/'run.json');d['source_extraction_sha256']=sha(p)
        failed_electric[name]=portable(d)
    if any(d['fit_accepted'] for d in failed_electric.values()):raise ValueError('Expected rejected electric attempts')
    weak=read(runs/'jtwpa4wm-circuit-refined/response.json')
    if weak['passed']:raise ValueError('Expected the under-resolved weak-idler failure')
    weak.pop('native_runs',None)
    model=load_circuit(kit/'backannotation.json');kr=read(kit/'backannotation.json')
    if kr['source']['gds_sha256']!=sha(kit/'kit.gds'):raise ValueError('KIT source changed')
    from twpa.models import TLS
    grid=np.linspace(0,8e9,801);linear=model.linear(grid)
    response=model.coupled([2e9],TLS(loss_tangent=1e-3,temperature=.02),frequency=grid,pump_frequency=6e9)
    thread=read(runs/'thread-control.json')
    if not thread['passed']:raise ValueError('Dense thread check failed')
    root=Path(__file__).resolve().parents[1]
    files=['src/twpa/em/'+p+'.py' for p in ('network','layout_annotation','backannotation','magnetic','magnetic_annotation','geometry','openems_worker','fasthenry','kit_annotation')]
    data=dict(schema_version=1,date=str(date.today()),code_version='0.10.0',
        source_file_sha256={p:sha(root/p) for p in files},cases=cases,
        kit=dict(backannotation=portable(kr),passive_max_energy_error=float(max(abs(abs(linear.s11)**2+abs(linear.s21)**2-1))),
                 coupled_gain=float(response.gain[0]),coupled_added_noise=float(response.added_noise[0]),pump_frequency_Hz=6e9,signal_frequency_Hz=2e9,
                 tls_assumptions=dict(loss_tangent=1e-3,temperature_K=.02,field_per_amp=0.)),
        rejected_coarse_magnetic=rejected,rejected_compact_snail=snail_rejected,rejected_electric=failed_electric,
        under_resolved_weak_idler=portable(weak),dense_thread_control=thread,
        scope='Hybrid individual PEC capacitor coupons and full finite-metal London MQS circuits; illustrative process. '
              'One-cell traveling-wave fixtures test finite assembly and conversion, not useful high gain. '
              'KIT is an explicit process/calibration periodic reduction, not full-layout EM. '
              'Classical extracted networks do not certify quantum noise or dynamical pump stability. '
              'Portable report only; source hashes do not transfer certificates to another GDS.')
    data['artifact_sha256']=digest(data)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Captured',len(cases),'accepted families:',output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--kit',type=Path,required=True);p.add_argument('--output',type=Path,default=Path('validation/layout-network-reference.json'))
    a=p.parse_args();capture(a.runs,a.kit,a.output)
