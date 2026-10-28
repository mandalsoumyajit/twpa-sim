"""Acceptance and back-annotation of full signed London current-loop matrices."""
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np

from ..layout.generate import digest
from ..layout.process import _json
from .backannotation import read_extraction
from .coupled_jpa import CoupledJPA,GalvanicJPA


def read_magnetic(value):
    d=deepcopy(value) if isinstance(value,dict) else _json(Path(value).read_text(encoding='utf-8'))
    h=d.pop('artifact_sha256',None)
    if h!=digest(d) or d.get('kind') not in ('jpa_london_network','layout_london_network'):raise ValueError('Magnetic artifact hash/type mismatch')
    d['artifact_sha256']=h
    return d


def _matrix(d):
    matrices=np.array([s['L_H'] for s in d['samples']])
    if not np.allclose(matrices,matrices[0],rtol=1e-7,atol=1e-20):raise ValueError('London lossless matrix is unexpectedly frequency dependent')
    return matrices[0]


def certify_magnetic(primary,controls):
    """Require independent fixed-terminal numerical controls.

    General terminal networks accept four total-L controls (bulk, film,
    contact, frequency). Separate energy decomposition additionally requires
    the six PEC/derivative controls used by the dedicated JPA workflow.

    London kinetic-energy matrix is dL/ds at s=1, with lambda -> sqrt(s)*lambda.
    For the stationary lossless quadratic current problem this is the kinetic
    energy bilinear form (Hellmann-Feynman). SC-minus-PEC is separately labeled.
    """
    p=read_magnetic(primary);l=_matrix(p);records={};loaded={}
    required={'bulk','thickness','contact','frequency','pec','pec_refined','plus','minus','plus_refined','minus_refined'}
    total_only=p['kind']=='layout_london_network' and set(controls)=={'bulk','thickness','contact','frequency'}
    if total_only:required=set(controls)
    if set(controls)!=required:raise ValueError('All declared magnetic convergence and kinetic-energy controls required')
    for name,v in controls.items():
        d=read_magnetic(v)
        for key in ('kind','source','manifest_sha256','stack','patches','reduction'):
            if d[key]!=p[key]:raise ValueError(f'Control changes physical {key}')
        # Mesh node numbers change on refinement. Compare physical finite
        # terminal identities/orientations, never their arbitrary node labels.
        legacy={'signal_path':('CAP_M2','JJ_M1'),'return_path':('JJ_M2','CAP_M1'),'pickup':('coil_lo','coil_hi')}
        def terminals(v):
            return [(q['name'],q['positive_patch'],q['negative_patch']) if 'positive_patch' in q else
                    (q['name'],*legacy[q['name']]) for q in v['ports']]
        if terminals(d)!=terminals(p):raise ValueError('Control changes physical port orientation')
        if d['run']['simulator']['binary_sha256']!=p['run']['simulator']['binary_sha256']:raise ValueError('Solver binary differs')
        loaded[name]=d
        if name!='frequency' and d['samples'][0]['frequency_Hz']!=p['samples'][0]['frequency_Hz']:raise ValueError('Controls must share the reference frequency')
    expected={'bulk':{'pitch_um'},'thickness':{'nh'},'contact':{'local_pitch_um'}}
    for name,keys in expected.items():
        d=loaded[name];changed={k for k in p['settings'] if k!='frequencies_hz' and d['settings'][k]!=p['settings'][k]}
        if not changed or not changed<=keys:raise ValueError('Control must change one numerical group')
        other=_matrix(d);delta=other-l
        diag=float(max(abs(np.diag(delta)/np.diag(l))))
        mutual=float(abs(delta[0,1])/max(abs(l[0,1]),1e-15)) if l.shape[0]>1 else 0
        if p['kind']=='layout_london_network':
            scale=np.sqrt(np.outer(np.diag(l),np.diag(l)));mask=~np.eye(len(l),dtype=bool)
            strong=mask & (abs(l)/scale>=.01)
            mutual=float(np.max(abs(delta[strong]/l[strong]))) if np.any(strong) else 0.
            coupling_error=float(np.max(abs(delta[mask])/scale[mask])) if np.any(mask) else 0.
            records[name]=dict(max_self_relative=diag,mutual_relative=mutual,coupling_absolute=coupling_error,
                               passed=diag<.02 and mutual<.05 and coupling_error<.005)
        else:records[name]=dict(max_self_relative=diag,mutual_relative=mutual,passed=diag<.02 and mutual<.05)
    if (p['settings']['pitch_um']>loaded['bulk']['settings']['pitch_um'] or
        p['settings']['nh']<loaded['thickness']['settings']['nh'] or
        p['settings']['local_pitch_um']>loaded['contact']['settings']['local_pitch_um']):raise ValueError('Primary must use finest declared numerical controls')
    fcontrol=loaded['frequency']
    if any(fcontrol['settings'][k]!=p['settings'][k] for k in p['settings'] if k!='frequencies_hz') or fcontrol['samples'][0]['frequency_Hz']<1e11:
        raise ValueError('Independent high-frequency control at fixed physical/current basis required')
    ferror=float(np.max(abs(_matrix(fcontrol)-l))/np.max(abs(l)))
    records['frequency']=dict(matrix_relative=ferror,passed=ferror<1e-7)
    if total_only:
        result=deepcopy(p);result.pop('artifact_sha256')
        result['primary_artifact_sha256']=p['artifact_sha256']
        result['certification']=dict(passed=all(r['passed'] for r in records.values()),checks=records,
            comparison_artifacts={k:d['artifact_sha256'] for k,d in loaded.items()},self_gate=.02,
            strong_mutual_gate=.05,strong_coupling_threshold=.01,absolute_coupling_gate=.005,
            scope='Total London network; no separate kinetic-energy or PEC-limit decomposition claimed')
        result['inductance_H']=l.tolist();result['artifact_sha256']=digest(result)
        return result
    for name in required-set(expected)-{'frequency'}:
        changed={k for k in p['settings'] if k!='frequencies_hz' and loaded[name]['settings'][k]!=p['settings'][k]}
        if changed!={'lambda_scale'}:raise ValueError('London controls must keep geometry/current mesh fixed')
    if p['settings']['lambda_scale']!=1:raise ValueError('Primary must use the declared operating-temperature lambdas')
    pec=_matrix(loaded['pec_refined']);pec_old=_matrix(loaded['pec'])
    if not 0<loaded['pec_refined']['settings']['lambda_scale']<loaded['pec']['settings']['lambda_scale']<=1e-3:raise ValueError('Invalid converged PEC limit')
    err=float(np.max(abs(pec-pec_old))/np.max(abs(pec)))
    records['pec_limit']=dict(matrix_relative=err,passed=err<1e-5,
                             model='Converged positive-lambda limit; exact lambda=0,sigma=0 is not a valid native conductor')
    derivatives=[]
    for pos,neg in [('plus','minus'),('plus_refined','minus_refined')]:
        sp=loaded[pos]['settings']['lambda_scale']**2;sm=loaded[neg]['settings']['lambda_scale']**2
        if not np.isclose((sp+sm)/2,1) or not 0<sp-sm<.5:raise ValueError('Centered lambda-squared derivative required')
        derivatives.append((_matrix(loaded[pos])-_matrix(loaded[neg]))/(sp-sm))
    kin=derivatives[1];derr=float(np.max(abs(derivatives[1]-derivatives[0]))/np.max(abs(kin)))
    mineig=float(np.linalg.eigvalsh((kin+kin.T)/2)[0])
    records['kinetic_derivative']=dict(matrix_relative=derr,minimum_eigenvalue_H=mineig,passed=derr<.02 and mineig>-1e-18)
    result=deepcopy(p);result.pop('artifact_sha256')
    result['primary_artifact_sha256']=p['artifact_sha256']
    result['topology']=p.get('topology','pickup' if len(p['ports'])==3 else 'galvanic')
    result['certification']=dict(passed=all(r['passed'] for r in records.values()),checks=records,
                                comparison_artifacts={k:d['artifact_sha256'] for k,d in loaded.items()},
                                self_gate=.02,mutual_gate=.05,derivative_gate=.02)
    result['inductance_H']=l.tolist();result['kinetic_energy_H']=kin.tolist()
    result['magnetic_energy_H']=(l-kin).tolist();result['pec_limit_H']=pec.tolist()
    result['SC_minus_PEC_H']=(l-pec).tolist()
    result['energy_note']='Only total London L is used in circuits. Kinetic energy is a converged lambda-squared variational derivative. SC-minus-PEC also includes magnetic redistribution.'
    result['artifact_sha256']=digest(result)
    return result


def annotate_coupled_jpa(gds,electric,magnetic,*,port_impedance=None,internal_resistance=1e8,junction_resistance=1e9,output=None):
    e=read_extraction(electric);m=read_magnetic(magnetic);path=Path(gds)
    if not e.get('convergence',{}).get('passed') or not m.get('certification',{}).get('passed'):raise ValueError('Certified electric and magnetic artifacts required')
    manifest=_json(path.with_suffix('.manifest.json').read_text(encoding='utf-8'))
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    for d in (e,m):
        if d['source']['gds_sha256']!=actual:raise ValueError('Stale extraction GDS')
        for key in ('design_sha256','process_sha256'):
            if d['source'][key]!=manifest[key]:raise ValueError('Stale process/design')
    mh=hashlib.sha256(path.with_suffix('.manifest.json').read_bytes()).hexdigest()
    if mh!=m['manifest_sha256'] or mh!=e['source']['manifest_sha256']:raise ValueError('Stale layout manifest')
    if not np.isclose(m['stack']['dielectric_gap_um'],e['stack']['dielectric_thickness_um'],rtol=0,atol=1e-12):raise ValueError('Electric/magnetic physical dielectric gap differs')
    l=np.array(m['inductance_H'])
    if l.shape!=(2,2):raise ValueError('A complete 2x2 extracted current matrix is required')
    pickup=manifest['design'].get('inductive_coupler',False)
    expected_topology='pickup' if pickup else 'galvanic'
    if m.get('topology',expected_topology)!=expected_topology:raise ValueError('Inductive reference-plane topology differs from layout')
    jj=e['source']['junction_estimates']
    cls=CoupledJPA if pickup else GalvanicJPA
    impedance=(50. if pickup else 5000.) if port_impedance is None else port_impedance
    c=cls(jj['Ic_A'],e['capacitance_F'],jj['Cj_F'],l[0,0],l[1,1],l[0,1],impedance,internal_resistance,junction_resistance)
    record=dict(schema_version=1,circuit_type=cls.__name__,circuit=asdict(c),topology=expected_topology,
                electric_sha256=e['artifact_sha256'],magnetic_sha256=m['artifact_sha256'],
                source=m['source'],manifest_sha256=mh,mode_parameters=c.modes(),
                certification=dict(electric=e['convergence'],magnetic=m['certification']),
                parameter_origins=dict(interconnect_inductance='Total London FastHenry closed tank-current loop, explicit return/vias' if pickup else 'Total London FastHenry JJ feed path and explicit common return',
                    pickup_inductance='Total London FastHenry physical pickup loop' if pickup else 'Total London FastHenry capacitor feed path and explicit common return',mutual_inductance='Signed native loop mutual; no fitted external kappa',
                    kinetic_inductance='Already included in total London matrix; energy decomposition is diagnostic only',
                    shunt_capacitance='openEMS actual capacitor electrodes',junction='Actual junction-window area times explicit process Ic/Cj density',
                    resistances='Explicit external embedding and user loss assumptions'),
                scope='Hybrid electric-capacitor and full-magnetic interconnect/pickup MQS network. No additional kinetic or empirical inductance is added.',
                limitations=['Local finite electrode feed fields excluded; explicit port-size sensitivity required for transfer to other partitions.',
                    'No dielectric loss, nonlinear kinetic response, radiation, or vortex/trapped-flux model.',
                    'Electric stray/cross capacitance outside the extracted MIM capacitor is not included; coupling here is the signed inductive coupling.'])
    if not pickup:
        origins=record['parameter_origins']
        origins['junction_feed_inductance']=origins.pop('interconnect_inductance')
        origins['capacitor_feed_inductance']=origins.pop('pickup_inductance')
    record['artifact_sha256']=digest(record)
    if output:Path(output).write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    return c,record
