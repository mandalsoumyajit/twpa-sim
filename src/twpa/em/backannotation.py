"""Validate extraction provenance, certify refinements, and annotate circuits."""
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from ..layout.generate import digest
from ..layout.process import _json, positive
from ..jpa_circuit import ParallelJPA


def read_extraction(path_or_dict):
    data=deepcopy(path_or_dict) if isinstance(path_or_dict,dict) else _json(Path(path_or_dict).read_text(encoding='utf-8'))
    recorded=data.pop('artifact_sha256',None)
    if recorded!=digest(data):
        raise ValueError('Extraction artifact hash mismatch')
    data['artifact_sha256']=recorded
    if data.get('schema_version')!=1 or data.get('kind') not in ('jpa_mim_capacitor','layout_mim_capacitor') or data.get('simulator',{}).get('name')!='openEMS':
        raise ValueError('Unsupported extraction artifact')
    if not data.get('fit_accepted'):
        raise ValueError('Extraction failed its fit/tail/passivity gates')
    positive(data['capacitance_F'],'extracted capacitance')
    if 'reused_native_extraction' in data:
        parent=data['reused_native_extraction']
        if 'reused_native_extraction' in parent:raise ValueError('Nested capacitor reuse is unsupported')
        parent=read_extraction(parent)
        for key in ('geometry','stack','stack_sha256','settings','simulator','capacitance_F','capacitance_by_frequency_F','convergence'):
            if data[key]!=parent[key]:raise ValueError('Reused native coupon changes physical/numerical evidence')
    return data


def reuse_capacitor(gds,component_id,certified):
    """Reuse a native result only for identical translated electrode geometry.

    The full original sealed artifact is retained. This is not a new field run,
    and a half-height end capacitor is never inferred from half an interior C.
    """
    from .geometry import capacitor_job
    from .stack import EMStack
    parent=read_extraction(certified)
    if 'reused_native_extraction' in parent:parent=parent['reused_native_extraction']
    if not isinstance(parent.get('convergence'),dict) or not parent['convergence'].get('passed'):
        raise ValueError('Certified native capacitor required')
    keys=('mesh_um','margin_um','z_cells','frequencies_hz','timesteps','port_resistance')
    settings={k:parent['settings'][k] for k in keys};settings['crop_margin_um']=parent['geometry']['crop_margin_um']
    job=capacitor_job(gds,EMStack(**parent['stack']),component_id=component_id,**settings)
    if job['geometry']!=parent['geometry'] or job['settings']!=parent['settings']:raise ValueError('Actual coupon differs; new electric extraction required')
    d=deepcopy(parent);d.pop('artifact_sha256');d.update(kind=job['kind'],source=job['source'],reused_native_extraction=parent)
    d['artifact_sha256']=digest(d)
    return d


def certify(primary, refinements, *, max_relative_change=.02):
    """Compare mesh, boundary-domain and duration sensitivity before annotation.

    Each refinement changes exactly one group of controls relative to primary:
    mesh=(mesh_um,z_cells), domain=margin_um, duration=timesteps. Comparison of
    coarse/fine spatial meshes is symmetric; labels retain both resolutions.
    An optional port_resistance comparison checks source-impedance invariance.
    Neither port-tail decay nor an RLC fit alone establishes convergence.
    """
    positive(max_relative_change,'max_relative_change')
    if max_relative_change>.05:
        raise ValueError('Convergence gate must be <=5%')
    p=read_extraction(primary); rows=[]; covered=set()
    groups={'mesh':{'mesh_um','z_cells'},'domain':{'margin_um'},'duration':{'timesteps'},'port_impedance':{'port_resistance'}}
    for value in refinements:
        r=read_extraction(value)
        for key in ('source','geometry','stack','stack_sha256','simulator'):
            if r[key]!=p[key]:
                raise ValueError(f'Convergence comparison changes physical {key}')
        if r['run']['worker_sha256']!=p['run']['worker_sha256']:
            raise ValueError('Convergence requires the same solver worker implementation')
        changed={key for key in p['settings'] if p['settings'][key]!=r['settings'][key]}
        names=[name for name,keys in groups.items() if changed and changed<=keys]
        if len(names)!=1:
            raise ValueError('Each convergence case must change only mesh, domain, duration or port impedance')
        name=names[0]
        a=np.asarray(p['capacitance_by_frequency_F']); b=np.asarray(r['capacitance_by_frequency_F'])
        change=float(np.max(abs(b/a-1)))
        rows.append(dict(control=name,relative_change=change,
                         primary_settings={k:p['settings'][k] for k in groups[name]},
                         comparison_settings={k:r['settings'][k] for k in groups[name]},
                         comparison_artifact_sha256=r['artifact_sha256'],passed=change<=max_relative_change))
        covered.add(name)
    if not {'mesh','domain','duration'}<=covered:
        raise ValueError('Mesh, domain and duration comparisons are all required')
    # Verify the final primary uses the finest/longest/largest compared controls.
    for value in refinements:
        r=read_extraction(value)
        if (p['settings']['mesh_um']>r['settings']['mesh_um'] or
            p['settings']['z_cells']<r['settings']['z_cells'] or
            p['settings']['margin_um']<r['settings']['margin_um'] or
            p['settings']['timesteps']<r['settings']['timesteps']):
            raise ValueError('Primary must use the finest mesh, largest domain and longest duration')
    result=deepcopy(p); result.pop('artifact_sha256')
    result['convergence']=dict(passed=all(r['passed'] for r in rows),
                               max_relative_change=max_relative_change,comparisons=rows)
    result['artifact_sha256']=digest(result)
    return result


def back_annotate_jpa(gds, extraction, *, port_impedance=5000., internal_resistance=1e7,
                      junction_resistance=1e9, output=None):
    """Replace the nominal shunt capacitor with a certified EM-extracted value.

    JJ area is measured from GDS; Ic/Cj densities and loss resistances remain
    explicit process/user assumptions. These are component-level circuits:
    interconnect/package/bias parasitics are not silently inferred from PEC.
    """
    data=read_extraction(extraction)
    convergence=data['convergence']
    if not isinstance(convergence,dict) or not convergence.get('passed'):
        raise ValueError('Mesh/domain/duration convergence certificate required')
    gds=Path(gds)
    manifest_path=gds.with_suffix('.manifest.json')
    m=_json(manifest_path.read_text(encoding='utf-8-sig'))
    if m['family']!='jpa':
        raise ValueError('Back-annotation supports single-junction jpa only')
    expected=dict(gds_sha256=hashlib.sha256(gds.read_bytes()).hexdigest(),
                  manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                  process_sha256=m['process_sha256'],design_sha256=m['design_sha256'])
    for key,value in expected.items():
        if value!=data['source'][key]:
            raise ValueError(f'Stale extraction: {key} differs from layout')
    jj=data['source']['junction_estimates']
    circuit=ParallelJPA(critical_current=jj['Ic_A'],shunt_capacitance=data['capacitance_F'],
                        junction_capacitance=jj['Cj_F'],junction_resistance=junction_resistance,
                        port_impedance=port_impedance,internal_resistance=internal_resistance)
    model=circuit.reduced_model()
    band=data['frequency_hz']
    if not band[0]<model.frequency<band[-1]:
        raise ValueError('Back-annotated resonance is outside the extraction fit band')
    record=dict(schema_version=1,circuit_type='ParallelJPA',circuit=asdict(circuit),
                reduced_model=asdict(model),source=expected,extraction_sha256=data['artifact_sha256'],
                stack_sha256=data['stack_sha256'],convergence=convergence,
                nominal_shunt_capacitance_F=data['source']['nominal_capacitance_F'],
                extracted_shunt_capacitance_F=data['capacitance_F'],
                parameter_origins=dict(shunt_capacitance='openEMS PEC capacitor coupon from GDS electrodes',
                                       critical_current='Actual GDS junction-window area times illustrative process Jc',
                                       junction_capacitance='Actual GDS junction-window area times illustrative process Cj density; added once',
                                       port_impedance='User-supplied external embedding',
                                       internal_resistance='User-supplied lumped loss; not extracted from PEC',
                                       junction_resistance='User-supplied RCSJ loss; not extracted from PEC'),
                limitations=['Capacitor coupon extraction, not a full-chip electromagnetic model.',
                             'Lead/loop/package inductance, kinetic inductance, cross-coupling, bias and port geometry remain unextracted.',
                             'Physical stack and junction densities are illustrative, not calibrated foundry data.'])
    record['artifact_sha256']=digest(record)
    if output:
        path=Path(output); path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return circuit,record


def load_circuit(path):
    """Load a sealed back-annotation artifact into its recorded circuit type."""
    data=_json(Path(path).read_text(encoding='utf-8'))
    stored=data.pop('artifact_sha256',None)
    if stored!=digest(data):
        raise ValueError('Back-annotation artifact hash mismatch')
    if data.get('schema_version')!=1:
        raise ValueError('Unsupported back-annotated circuit')
    if data.get('circuit_type')=='KITPeriodicModel':
        from .kit_annotation import load_kit_data
        if not data.get('checks',{}).get('geometry_passed'):raise ValueError('Checked KIT geometry required')
        return load_kit_data(data['circuit'])
    if data.get('circuit_type')=='LayoutNetwork':
        from .network import LayoutNetwork
        c=data.get('certification',{})
        if not c.get('magnetic',{}).get('passed') or not c.get('electric') or not all(v.get('passed') for v in c['electric'].values()):
            raise ValueError('Complete certified extraction required')
        return LayoutNetwork(data['circuit'])
    if data.get('circuit_type') in ('CoupledJPA','GalvanicJPA'):
        from .coupled_jpa import CoupledJPA,GalvanicJPA
        if not all(data.get('certification',{}).get(k,{}).get('passed') for k in ('electric','magnetic')):
            raise ValueError('Certified electric and magnetic extraction required')
        cls=CoupledJPA if data['circuit_type']=='CoupledJPA' else GalvanicJPA
        circuit=cls(**data['circuit'])
        if not np.allclose([m['frequency'] for m in circuit.modes()],[m['frequency'] for m in data['mode_parameters']],rtol=1e-12):
            raise ValueError('Stored circuit modes differ')
        return circuit
    if data.get('circuit_type')!='ParallelJPA':raise ValueError('Unsupported back-annotated circuit')
    if not data.get('convergence',{}).get('passed'):
        raise ValueError('Uncertified extraction')
    circuit=ParallelJPA(**data['circuit'])
    if not np.isclose(circuit.reduced_model().frequency,data['reduced_model']['frequency'],rtol=1e-12):
        raise ValueError('Stored reduced model differs from circuit')
    return circuit
