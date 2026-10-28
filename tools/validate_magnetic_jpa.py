"""Complete London-interconnect/pickup JPA extraction and JoSIM comparison.

Fresh run directories only. FastHenry uses the selected installed Whiteley
binary with native wall/memory limits. Raw JoSIM waveforms are discarded
after full reads; circuit inputs, logs, hashes and compact results remain.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np

from twpa.layout import load_process
from twpa.layout.generate import generate,JosephsonLayout
from twpa.em import OpenEMS,load_stack,capacitor_job,certify
from twpa.em.fasthenry import FastHenry
from twpa.em.magnetic import jpa_magnetic_job,load_magnetic_stack
from twpa.em.magnetic_annotation import read_magnetic,certify_magnetic,annotate_coupled_jpa
from twpa.em.coupled_jpa import pump,response


def extract(gds,root,runner,*,primary_artifact=None,workers=2):
    stack=load_magnetic_stack()
    if workers not in (1,2):raise ValueError('Only one or two bounded extraction workers are supported')
    base=dict(pitch_um=1.5,local_pitch_um=.25,nh=3,frequencies_hz=(1e6,))
    if primary_artifact:
        primary=read_magnetic(primary_artifact)
        job=jpa_magnetic_job(gds,stack,**base)
        if primary['source']!=job['source'] or primary['stack']!=job['stack']:raise ValueError('Provided primary is stale')
        if primary['settings']!=job['settings']:raise ValueError('Provided primary does not match production numerical controls')
    else:
        job=jpa_magnetic_job(gds,stack,**base)
        print('FastHenry primary',job['audit']['filaments'],flush=True)
        primary=runner.run(job,root/'magnetic-primary')
    cases=dict(bulk=dict(pitch_um=2.),thickness=dict(nh=1),contact=dict(local_pitch_um=.5),frequency=dict(frequencies_hz=(1e11,)),
               pec=dict(lambda_scale=1e-3),pec_refined=dict(lambda_scale=1e-4),
               plus=dict(lambda_scale=float(np.sqrt(1.1))),minus=dict(lambda_scale=float(np.sqrt(.9))),
               plus_refined=dict(lambda_scale=float(np.sqrt(1.05))),minus_refined=dict(lambda_scale=float(np.sqrt(.95))))
    results={}
    def one(name,settings):
        job=jpa_magnetic_job(gds,stack,**(base|settings))
        print('FastHenry',name,job['audit']['filaments'],flush=True)
        value=runner.run(job,root/('magnetic-'+name))
        print(name,'L (pH)',np.array(value['samples'][0]['L_H'])*1e12,flush=True)
        return value
    with ThreadPoolExecutor(max_workers=workers) as pool:
        pending={pool.submit(one,name,settings):name for name,settings in cases.items()}
        for future in as_completed(pending):results[pending[future]]=future.result()
    certificate=certify_magnetic(primary,results)
    (root/'magnetic-certified.json').write_text(json.dumps(certificate,indent=2)+'\n',encoding='utf-8')
    if not certificate['certification']['passed']:raise RuntimeError('Magnetic acceptance failed; inspect all retained controls')
    return certificate


def extract_cap(gds,root,runner):
    cases=dict(primary=dict(mesh_um=.5,z_cells=3,margin_um=7,timesteps=100000),
               mesh=dict(mesh_um=1,z_cells=2,margin_um=7,timesteps=100000),
               domain=dict(mesh_um=.5,z_cells=3,margin_um=5,timesteps=100000),
               duration=dict(mesh_um=.5,z_cells=3,margin_um=7,timesteps=50000))
    results={}
    for name,settings in cases.items():
        print('openEMS',name,flush=True)
        results[name]=runner.run(capacitor_job(gds,load_stack(),**settings),root/('electric-'+name))
    d=certify(results['primary'],[results[n] for n in ('mesh','domain','duration')])
    (root/'electric-certified.json').write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8')
    if not d['convergence']['passed']:raise RuntimeError('Electric refinement gates failed')
    return d


def select_pump(circuit):
    mode=circuit.modes()[0];f=mode['frequency'];gamma=mode['kappa']/(2*np.pi)
    for detuning in [-.5,-1.,0.]:
        fp=f+detuning*gamma;fs=fp+gamma*np.array([.25,.5]);previous=None
        for power in np.geomspace(2e-16,2e-12,33):
            try:
                point=pump(circuit,fp,float(power),harmonics=3,initial=previous)
                previous=point;result=response(point,fs)
            except (ValueError,RuntimeError):continue
            if max(result['gain'])>=2. and point.phase_peak<.6:return point,fs,result
    raise RuntimeError('No stable weakly nonlinear amplifier point in declared search range')


def compare(circuit,root,executable,wsl=None,corrector=0,step=5e-15,start=700e-9,stop=1300e-9,trace_format='npy'):
    from twpa.backends import JoSIM
    from twpa.jpa_circuit import Tone,Transient
    from twpa.josim_jpa import run_circuit,sideband_measurement,waveforms
    from twpa.measurement import fit_tones
    runner=JoSIM(executable,wsl,timeout=240,corrector=corrector,output_format=trace_format)
    transient=Transient(step,stop,start,2e-13 if trace_format=='npy' else 1e-12);mode=circuit.modes()[0]
    f=mode['frequency'];gamma=mode['kappa']/(2*np.pi)
    point,fs,predicted=select_pump(circuit)
    fitted=np.unique(np.r_[fs,2*point.frequency-fs,point.frequency,2*point.frequency,3*point.frequency])
    required=1.1/np.min(np.diff(fitted))
    if stop-start<required:raise ValueError(f'Acquisition requires at least {required:g} s to resolve the selected pump/probe/idler tones')
    def run(label,tones):
        result=run_circuit(runner,circuit,root/label,tones,transient)
        result.waveform_path.unlink()
        return result
    passive_fs=f+gamma*np.array([-.5,0,.5]);probe=1e-23
    actual=run('passive-probe',[Tone(v,probe) for v in passive_fs])
    incoming,outgoing=waveforms(actual,circuit.port_impedance)
    s=fit_tones(actual.time,outgoing,passive_fs)/fit_tones(actual.time,incoming,passive_fs)
    measured=dict(signal_frequency=passive_fs.tolist(),signal_real=s.real.tolist(),signal_imag=s.imag.tolist(),gain=(abs(s)**2).tolist())
    error=float(max(abs(s-circuit.passive(passive_fs))))
    base=run('pump-base',[Tone(point.frequency,point.power)])
    actual=run('pump-probe',[Tone(point.frequency,point.power)]+[Tone(v,probe) for v in fs])
    pm=sideband_measurement(actual,base,circuit.port_impedance,point.frequency,fs)
    gain_error=10*np.log10(np.array(pm['gain'])/predicted['gain'])
    idler_error=10*np.log10(np.array(pm['idler_conversion'])/predicted['idler_conversion'])
    refined=pump(circuit,point.frequency,point.power,harmonics=4,initial=point)
    high=response(refined,fs,sideband_order=3)
    hc_gain=float(max(abs(10*np.log10(high['gain']/predicted['gain']))))
    hc_idler=float(max(abs(10*np.log10(high['idler_conversion']/predicted['idler_conversion']))))
    passed=bool(error<.02 and max(abs(gain_error))<.25 and max(abs(idler_error))<.25 and hc_gain<.02 and hc_idler<.02)
    result=dict(passed=passed,simulator=runner.probe(),circuit=asdict(circuit),modes=circuit.modes(),
                transient=asdict(transient),passive=dict(measured=measured,max_complex_error=error,
                    model_real=circuit.passive(passive_fs).real.tolist(),model_imag=circuit.passive(passive_fs).imag.tolist()),
                pump=dict(frequency=point.frequency,power=point.power,phase_peak=point.phase_peak,
                          floquet_radius=point.floquet_radius,residual=point.residual),
                pumped=dict(measured=pm,model_gain=predicted['gain'].tolist(),model_idler_conversion=predicted['idler_conversion'].tolist(),
                            gain_error_db=gain_error.tolist(),idler_error_db=idler_error.tolist()),
                harmonic_refinement=dict(gain_change_db=hc_gain,idler_change_db=hc_idler),
                acceptance=dict(complex_error=.02,gain_error_db=.25,idler_error_db=.25,harmonic_change_db=.02))
    if wsl:
        binary=subprocess.run(['wsl','-d',wsl,'--exec','which',executable],capture_output=True,text=True,check=True).stdout.strip()
        result['simulator_binary_sha256']=subprocess.run(['wsl','-d',wsl,'--exec','sha256sum',binary],capture_output=True,text=True,check=True).stdout.split()[0]
    else:result['simulator_binary_sha256']=hashlib.sha256(Path(shutil.which(executable) or executable).read_bytes()).hexdigest()
    (root/'response.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result['pumped']),flush=True)
    if not passed:raise RuntimeError('Coupled amplifier/circuit comparison failed; refine or diagnose before acceptance')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--gds',type=Path)
    p.add_argument('--magnetic-primary',type=Path);p.add_argument('--installation',type=Path)
    p.add_argument('--workers',type=int,default=2,choices=(1,2),help='Each native case retains its 180s/4GiB cap; max aggregate address-space request 8GiB')
    p.add_argument('--openems-python');p.add_argument('--openems-dll-directory');p.add_argument('--extract-magnetic-only',action='store_true')
    p.add_argument('--electric',type=Path);p.add_argument('--magnetic',type=Path)
    p.add_argument('--executable',default='josim-cli');p.add_argument('--wsl');p.add_argument('--corrector',type=int,default=0)
    p.add_argument('--trace-format',choices=('npy','csv'),default='npy')
    p.add_argument('--step',type=float,default=5e-15);p.add_argument('--start',type=float,default=700e-9);p.add_argument('--stop',type=float,default=1300e-9)
    p.add_argument('--skip-josim',action='store_true')
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    gds=root/'jpa.gds'
    if a.gds:
        shutil.copy2(a.gds,gds);shutil.copy2(a.gds.with_suffix('.manifest.json'),gds.with_suffix('.manifest.json'))
    else:generate(JosephsonLayout(inductive_coupler=True),load_process('jj_overlap_illustrative'),gds)
    magnetic=read_magnetic(a.magnetic) if a.magnetic else extract(gds,root,FastHenry.from_installation(a.installation),primary_artifact=a.magnetic_primary,workers=a.workers)
    if a.extract_magnetic_only:return
    if a.electric:electric=json.loads(a.electric.read_text())
    else:
        if not a.openems_python:p.error('--openems-python or a matching certified --electric artifact required')
        electric=extract_cap(gds,root,OpenEMS(a.openems_python,a.openems_dll_directory))
    circuit,record=annotate_coupled_jpa(gds,electric,magnetic,output=root/'backannotation.json')
    if not a.skip_josim:compare(circuit,root,a.executable,a.wsl,a.corrector,a.step,a.start,a.stop,a.trace_format)


if __name__=='__main__':main()
