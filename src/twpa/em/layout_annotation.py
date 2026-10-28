"""Bind every extracted component and finite metal terminal to an actual GDS."""
import hashlib
import json
from pathlib import Path
import numpy as np
from ..layout.generate import digest,database
from ..layout.process import _json,Process
from .backannotation import read_extraction
from .magnetic_annotation import read_magnetic
from .magnetic import layout_magnetic_job,MagneticStack
from .network import LayoutNetwork


def back_annotate_layout(gds,electric,magnetic,*,port_impedance=50.,junction_resistance=1e4,
                         dielectric_resistance=1e4,output=None):
    """Return a nonuniform terminal circuit; require all component certificates.

    No uniform-cell fit, loop-L replacement or extra kinetic inductance is
    performed. IND annotations already belong to the extracted metal network.
    """
    path=Path(gds);mp=path.with_suffix('.manifest.json');m=_json(mp.read_text(encoding='utf-8'))
    h=hashlib.sha256(path.read_bytes()).hexdigest();mh=hashlib.sha256(mp.read_bytes()).hexdigest()
    field=read_magnetic(magnetic)
    if field['kind']!='layout_london_network' or not isinstance(field.get('certification'),dict) or not field['certification'].get('passed'):
        raise ValueError('Certified full terminal magnetic network required')
    def bound(d):
        s=d['source']
        if s['gds_sha256']!=h or s['process_sha256']!=m['process_sha256'] or s['design_sha256']!=m['design_sha256']:
            raise ValueError('Stale layout/process/design extraction')
        actual=d.get('manifest_sha256',s.get('manifest_sha256'))
        if actual!=mh:raise ValueError('Stale manifest extraction')
    bound(field)
    job=layout_magnetic_job(path,MagneticStack(**field['stack']),**field['settings'],max_filaments=max(20000,field['audit']['filaments']))
    primary_hash=field.get('primary_artifact_sha256',field['artifact_sha256'])
    if field['run']['job_sha256']!=digest(job):raise ValueError('Current geometry/terminal basis differs from extraction job')
    caps={c['id']:c for c in m['components'] if c['kind']=='CAP'}
    if not isinstance(electric,dict) or set(electric)!=set(caps):raise ValueError('One electric certificate per actual capacitor ID required')
    capacitor_records=[];hashes={}
    for name in caps:
        d=read_extraction(electric[name]);bound(d)
        if d['source']['component_id']!=name or not isinstance(d.get('convergence'),dict) or not d['convergence'].get('passed'):
            raise ValueError('Missing/wrong capacitor convergence certificate')
        if not np.isclose(d['stack']['dielectric_thickness_um'],field['stack']['dielectric_gap_um'],rtol=0,atol=1e-12):
            raise ValueError('Electric/magnetic physical gap differs')
        hashes[name]=d['artifact_sha256']
        capacitor_records.append(dict(id=name,positive=name+'_M2',negative=name+'_M1',
            capacitance=d['capacitance_F'],resistance=dielectric_resistance,
            extraction_frequency_band_Hz=[min(d['settings']['frequencies_hz']),max(d['settings']['frequencies_hz'])]))
    process=Process(m['process']);db=database();layout=db.Layout();layout.read(str(path));top=layout.top_cell()
    regions=db.Region(top.begin_shapes_rec(layout.layer(*process.layer('JJ'))))
    junctions=[]
    for c in m['components']:
        if c['kind']!='JJ':continue
        box=db.Box(*(process.coordinate(v) for v in c['box_um']))
        area=(regions&db.Region(box)).area()*layout.dbu**2
        if area<=0:raise ValueError('Missing actual junction window')
        junctions.append(dict(id=c['id'],positive=c['id']+'_M1',negative=c['id']+'_M2',
            area_um2=area,critical_current=area*process.parameter('critical_current_density','uA/um2')*1e-6,
            capacitance=area*process.parameter('junction_specific_capacitance','fF/um2')*1e-15,
            resistance=junction_resistance))
    ports=[dict(positive='INPUT',negative='GROUND0',impedance=port_impedance)]
    if 'OUTPUT' in {p['name'] for p in job['patches']}:ports.append(dict(positive='OUTPUT',negative='GROUND1',impedance=port_impedance))
    branches=[dict(positive=p['positive_patch'],negative=p['negative_patch']) for p in field['ports']]
    data=dict(family=m['family'],nodes=[p['name'] for p in job['patches']],reference='GROUND0',
              inductance_H=field['inductance_H'],metal_branches=branches,junctions=junctions,
              capacitors=capacitor_records,ports=ports)
    circuit=LayoutNetwork(data)
    record=dict(schema_version=1,circuit_type='LayoutNetwork',circuit=data,source=field['source'],manifest_sha256=mh,
                magnetic_sha256=field['artifact_sha256'],primary_magnetic_sha256=primary_hash,electric_sha256=hashes,
                certification=dict(magnetic=field['certification'],electric={k:read_extraction(v)['convergence'] for k,v in electric.items()}),
                loops=circuit.cycles.tolist(),
                parameter_origins=dict(metal='Full signed finite-conductor London matrix; kinetic energy already included',
                    capacitors='Individual certified openEMS coupons including distinct end geometry',
                    junctions='Actual window areas times explicit process densities; Cj added once',
                    losses='User-specified constant resistances; not extracted from PEC/lossless London',
                    flux='Externally prescribed signed flux per superconducting loop; no extracted bias-coil mask'),
                limitations=['Hybrid component electric and full finite-metal London MQS circuit, not full electric cross-capacitance extraction.',
                    'Explicit DC branch and AC phase guard; harmonic-balance convergence alone does not certify dynamical pump stability.',
                    'Illustrative material/density assumptions; no fabrication qualification or nonlinear kinetic/phase-slip model.'])
    record['artifact_sha256']=digest(record)
    if output:
        dest=Path(output);dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return circuit,record
