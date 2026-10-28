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
    if data.get('schema_version')!=1 or data.get('kind')!='jpa_mim_capacitor' or data.get('simulator',{}).get('name')!='openEMS':
        raise ValueError('Unsupported extraction artifact')
    if not data.get('fit_accepted'):
        raise ValueError('Extraction failed its fit/tail/passivity gates')
    positive(data['capacitance_F'],'extracted capacitance')
    return data


def certify(primary, refinements, *, max_relative_change=.02):
    """Compare mesh, boundary-domain and duration sensitivity before annotation.

    Each refinement changes exactly one group of controls relative to primary:
    mesh=(mesh_um,z_cells), domain=margin_um, duration=timesteps. Comparison of
    coarse/fine spatial meshes is symmetric; labels retain both resolutions.
    Neither port-tail decay nor an RLC fit alone establishes convergence.
    """
    positive(max_relative_change,'max_relative_change')
    if max_relative_change>.05:
        raise ValueError('Convergence gate must be <=5%')
    p=read_extraction(primary); rows=[]; covered=set()
    groups={'mesh':{'mesh_um','z_cells'},'domain':{'margin_um'},'duration':{'timesteps'}}
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
            raise ValueError('Each convergence case must change only mesh, domain or duration')
        name=names[0]
        a=np.asarray(p['capacitance_by_frequency_F']); b=np.asarray(r['capacitance_by_frequency_F'])
        change=float(np.max(abs(b/a-1)))
        rows.append(dict(control=name,relative_change=change,
                         primary_settings={k:p['settings'][k] for k in groups[name]},
                         comparison_settings={k:r['settings'][k] for k in groups[name]},
                         comparison_artifact_sha256=r['artifact_sha256'],passed=change<=max_relative_change))
        covered.add(name)
    if covered!=set(groups):
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
    """Load a checked back-annotation artifact for existing JPA/JoSIM APIs."""
    data=_json(Path(path).read_text(encoding='utf-8'))
    stored=data.pop('artifact_sha256',None)
    if stored!=digest(data):
        raise ValueError('Back-annotation artifact hash mismatch')
    if data.get('schema_version')!=1 or data.get('circuit_type')!='ParallelJPA':
        raise ValueError('Unsupported back-annotated circuit')
    if not data.get('convergence',{}).get('passed'):
        raise ValueError('Uncertified extraction')
    circuit=ParallelJPA(**data['circuit'])
    if not np.isclose(circuit.reduced_model().frequency,data['reduced_model']['frequency'],rtol=1e-12):
        raise ValueError('Stored reduced model differs from circuit')
    return circuit
