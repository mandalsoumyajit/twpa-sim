"""Extract individual layout components and compare a finite network with JoSIM.

Each run retains source/material snapshots and raw signed matrices. Only
geometrically identical translated capacitor coupons may reuse native evidence.
External flux is prescribed per physical loop; bias-coil geometry is not inferred.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
from twpa.layout import load_process
from twpa.layout.generate import generate,JosephsonLayout
from twpa.em.magnetic import layout_magnetic_job
from twpa.em.fasthenry import FastHenry
from twpa.em.magnetic_annotation import read_magnetic,certify_magnetic
from twpa.em import OpenEMS,load_stack,capacitor_job,certify
from twpa.em.backannotation import read_extraction,reuse_capacitor
from twpa.em.layout_annotation import back_annotate_layout


def save(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def magnetic(gds,folder,runner,primary=None,workers=1,*,pitch_um=2.,local_pitch_um=.125,patch_side_um=1.,nh=2,max_filaments=40000,contact_pitch_um=None):
    if type(nh) is not int or nh<2 or workers not in (1,2):raise ValueError('Use nh >=2 and one or two workers')
    if contact_pitch_um is not None and contact_pitch_um<=local_pitch_um:raise ValueError('Contact control must be coarser than primary')
    base=dict(pitch_um=pitch_um,local_pitch_um=local_pitch_um,patch_side_um=patch_side_um,nh=nh,frequencies_hz=(1e6,))
    if primary:
        p=read_magnetic(primary)
        if p['settings']!={**base,'frequencies_hz':[1e6],'lambda_scale':1.}:
            raise ValueError('Supplied primary has different numerical settings')
        if p['run']['job_sha256']!=__import__('twpa.layout.generate',fromlist=['digest']).digest(layout_magnetic_job(gds,**p['settings'],max_filaments=max_filaments)):
            raise ValueError('Supplied primary differs from current source/job')
    else:p=runner.run(layout_magnetic_job(gds,**base,max_filaments=max_filaments),folder/'magnetic-primary')
    cases=dict(bulk=dict(pitch_um=1.5*pitch_um),thickness=dict(nh=max(1,nh-1)),contact=dict(local_pitch_um=contact_pitch_um or 2*local_pitch_um),frequency=dict(frequencies_hz=(1e11,)))
    def one(name,overrides):
        print('FastHenry',name,flush=True)
        return runner.run(layout_magnetic_job(gds,**(base|overrides),max_filaments=max_filaments),folder/('magnetic-'+name))
    controls={}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        pending={pool.submit(one,n,s):n for n,s in cases.items()}
        for future in as_completed(pending):controls[pending[future]]=future.result()
    result=certify_magnetic(p,controls);save(folder/'magnetic-certified.json',result)
    if not result['certification']['passed']:raise RuntimeError('Magnetic convergence failed; retain gates and refine')
    return result


def electric(gds,folder,runner,library=(),*,port_resistance=1.,port_control_resistance=None,timesteps=100000):
    if type(timesteps) is not int or timesteps<20000:raise ValueError('At least 20,000 primary timesteps required')
    manifest=json.loads(gds.with_suffix('.manifest.json').read_text());cache=[read_extraction(p) for p in library];results={}
    settings=dict(primary=dict(mesh_um=.5,z_cells=3,margin_um=7,timesteps=100000),
        mesh=dict(mesh_um=1.,z_cells=2,margin_um=7,timesteps=100000),
        domain=dict(mesh_um=.5,z_cells=3,margin_um=5,timesteps=100000),
        duration=dict(mesh_um=.5,z_cells=3,margin_um=7,timesteps=50000))
    for setting in settings.values():
        setting['frequencies_hz']=(1e9,2e9,3e9,4e9,5e9,6e9,10e9,20e9,40e9,60e9,80e9,100e9)
        setting['port_resistance']=port_resistance
        setting['timesteps']=timesteps
    settings['duration']['timesteps']=timesteps//2
    if port_control_resistance is not None:
        if port_control_resistance==port_resistance:raise ValueError('Port-impedance control must differ')
        settings['port_impedance']=settings['primary']|dict(port_resistance=port_control_resistance)
    for component in manifest['components']:
        if component['kind']!='CAP':continue
        name=component['id'];job=capacitor_job(gds,load_stack(),component_id=name,**settings['primary'])
        candidates=[d for d in cache if d['geometry']==job['geometry'] and d['stack']==job['stack'] and d['settings']==job['settings']
                    and isinstance(d.get('convergence'),dict) and d['convergence'].get('passed')]
        if port_control_resistance is not None:
            candidates=[d for d in candidates if any(row.get('control')=='port_impedance' and
                row.get('comparison_settings',{}).get('port_resistance')==port_control_resistance
                for row in d.get('convergence',{}).get('comparisons',[]))]
        if candidates:
            results[name]=reuse_capacitor(gds,name,candidates[0]);print('Reuse identical native coupon',name,flush=True)
        else:
            if runner is None:raise ValueError('A new electrode footprint needs an openEMS runtime')
            cases={}
            for control,s in settings.items():
                print('openEMS',name,control,flush=True)
                cases[control]=runner.run(capacitor_job(gds,load_stack(),component_id=name,**s),folder/('electric-'+name+'-'+control))
            results[name]=certify(cases['primary'],[d for c,d in cases.items() if c!='primary'])
            if not results[name]['convergence']['passed']:raise RuntimeError('Electric convergence failed')
            cache.append(results[name])
        save(folder/('electric-'+name+'-certified.json'),results[name])
    return results


def measure(run,circuit,frequencies,extra=(),start=None):
    from twpa.measurement import fit_tones
    def voltage(name):return np.zeros_like(run.time) if name==circuit.reference else run.trace('V(T'+str(circuit.index[name])+')')
    source=run.trace('V(SRC0)')
    fitted=list(frequencies)+[float(f) for f in extra if f not in frequencies]
    mask=np.ones(len(run.time),bool) if start is None else run.time>=start
    def fit(value):return fit_tones(run.time[mask],value[mask],fitted)[:len(frequencies)]
    inc=source/(2*np.sqrt(circuit.R[0]));outputs=[]
    for i,p in enumerate(circuit.data['ports']):
        value=(voltage(p['positive'])-voltage(p['negative']))/np.sqrt(circuit.R[i])
        if i==0:value=value-inc
        outputs.append(fit(value))
    return fit(inc),np.array(outputs)


def compare(circuit,folder,runner,*,pump_frequency=5e9,signal_frequency=None,pump_power=1e-10,
            flux_bias=.05,flux_pump=0.,order=2,step=5e-15,start=80e-9,stop=160e-9):
    from twpa.jpa_circuit import Tone,Transient
    from twpa.measurement import fit_tones
    simulator=dict(runner.probe())
    if runner.wsl_distribution:
        simulator['binary_sha256']=subprocess.run(['wsl','-d',runner.wsl_distribution,'--exec','sha256sum',runner.executable],
            capture_output=True,text=True,check=True,timeout=30).stdout.split()[0]
    else:simulator['binary_sha256']=hashlib.sha256(Path(shutil.which(runner.executable) or runner.executable).read_bytes()).hexdigest()
    nl=len(circuit.cycles);bias=flux_bias if nl else 0.
    fs=signal_frequency or (.4*pump_frequency if order==1 else .8*pump_frequency);fi=order*pump_frequency-fs
    if fi<=0 or abs(fi-fs)<1:raise ValueError('Distinct positive signal and idler required')
    point=circuit.periodic(pump_frequency,pump_power,flux_bias=bias,flux_pump=flux_pump,harmonics=4)
    predicted=circuit.response(point,[fs],sidebands=3)
    refined=circuit.periodic(pump_frequency,pump_power,flux_bias=bias,flux_pump=flux_pump,harmonics=6,quadrature=512)
    high=circuit.response(refined,[fs],sidebands=4,quadrature=1024)
    target_port=len(circuit.R)-1;bands=predicted['bands'];s=predicted['waves'][0,bands==0,target_port][0]
    idler=np.conj(predicted['waves'][0,bands==-order,target_port][0])*np.sqrt(fs/fi)
    hs=high['waves'][0,high['bands']==0,target_port][0];hi=np.conj(high['waves'][0,high['bands']==-order,target_port][0])*np.sqrt(fs/fi)
    spectral=dict(signal_change_db=float(20*np.log10(abs(hs/s))),idler_change_db=float(20*np.log10(abs(hi/idler))))
    transient=Transient(step,stop,start,1e-12);probe=pump_power*1e-5 if pump_power else 1e-15
    native={}
    def run(name,tones,ac_flux=flux_pump):
        deck=circuit.netlist(tones,transient,flux_bias=bias,flux_pump=ac_flux,flux_frequency=pump_frequency)
        result=runner.run(deck,folder/name,metadata=dict(flux_bias=bias,flux_pump=ac_flux,circuit=circuit.data))
        native[name]=result.manifest
        return result
    passive_fs=np.array([fs,fi]);passive=run('passive',[Tone(float(f),probe) for f in passive_fs],0.)
    incoming,outgoing=measure(passive,circuit,passive_fs);actual_passive=(outgoing/incoming).T
    passive_error=float(np.max(abs(actual_passive-circuit.passive(passive_fs,flux_bias=bias))))
    tones=[Tone(pump_frequency,pump_power,ramp_time=10e-9)]
    base=run('pump-base',tones);actual=run('pump-probe',tones+[Tone(fs,probe)])
    if not np.array_equal(actual.time,base.time):raise ValueError('Pump baseline and probe timestamps differ')
    # Measure baseline-subtracted voltage directly, avoiding division by the
    # absent idler input. Incoming signal is the source wave at fs.
    nuisance=np.arange(1,7)*pump_frequency
    a,out=measure(actual,circuit,[fs,fi],nuisance);_,baseline=measure(base,circuit,[fs,fi],nuisance)
    # Fit also pump harmonics to suppress leakage when acquisition is not an
    # exact number of periods; the selected GHz probes are well resolved.
    measured_signal=(out-baseline)[target_port,0]/a[0]
    measured_idler=(out-baseline)[target_port,1]/np.conj(a[0])*np.sqrt(fs/fi)
    late=(start+stop)/2
    la,lo=measure(actual,circuit,[fs,fi],nuisance,start=late)
    _,lb=measure(base,circuit,[fs,fi],nuisance,start=late)
    late_signal=(lo-lb)[target_port,0]/la[0]
    late_idler=(lo-lb)[target_port,1]/np.conj(la[0])*np.sqrt(fs/fi)
    acquisition=dict(late_start=late,signal_change_db=float(20*np.log10(abs(late_signal/measured_signal))),
                     idler_change_db=float(20*np.log10(abs(late_idler/measured_idler))))
    gain_error=float(20*np.log10(abs(measured_signal/s)));idler_error=float(20*np.log10(abs(measured_idler/idler)))
    dc_model=circuit.BJ.T@point['dc'];dc_actual=np.array([np.mean(base.trace('P(BJJ'+str(i)+')')) for i in range(len(circuit.ic))])
    dc_error=float(np.max(abs(dc_actual-dc_model)))
    result=dict(passed=bool(passive_error<.02 and abs(gain_error)<.25 and abs(idler_error)<.25 and dc_error<1e-4 and
                          max(abs(v) for v in spectral.values())<.02 and
                          max(abs(acquisition[k]) for k in ('signal_change_db','idler_change_db'))<.02),
        source='Actual native JoSIM traces; model coefficients not fitted',simulator=simulator,native_runs=native,
        family=circuit.data['family'],circuit=circuit.data,
        drive=dict(pump_frequency=pump_frequency,signal_frequency=fs,idler_frequency=fi,pump_power=pump_power,
                   probe_power=probe,flux_bias=bias,flux_pump=flux_pump,mixing_order=order),
        transient=dict(step=step,start=start,stop=stop,output_step=1e-12),
        passive=dict(max_complex_error=passive_error,model_real=circuit.passive(passive_fs,flux_bias=bias).real.tolist(),
                     model_imag=circuit.passive(passive_fs,flux_bias=bias).imag.tolist(),measured_real=actual_passive.real.tolist(),measured_imag=actual_passive.imag.tolist()),
        pumped=dict(model_gain=abs(s)**2,measured_gain=abs(measured_signal)**2,
                    model_idler_conversion=abs(idler)**2,measured_idler_conversion=abs(measured_idler)**2,
                    gain_error_db=gain_error,idler_error_db=idler_error,maximum_phase=point['maximum_phase'],residual=point['residual'],
                    junction_dc_model=dc_model.tolist(),junction_dc_measured=dc_actual.tolist(),max_dc_phase_error=dc_error),
        spectral_refinement=spectral,acquisition_refinement=acquisition,
        acceptance=dict(passive_complex=.02,gain_db=.25,idler_db=.25,spectral_db=.02,acquisition_db=.02))
    save(folder/'response.json',result);print(json.dumps(result['pumped']),flush=True)
    if not result['passed']:raise RuntimeError('Circuit/JoSIM comparison failed; inspect and refine')
    for run_result in (passive,base,actual):run_result.waveform_path.unlink()
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',choices=['squid_jpa','jtwpa4wm','rfsquid3wm','snail3wm'],required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--gds',type=Path);p.add_argument('--cells',type=int,default=1)
    p.add_argument('--large-junctions',type=int,default=3)
    p.add_argument('--installation',type=Path);p.add_argument('--magnetic-primary',type=Path);p.add_argument('--magnetic',type=Path)
    p.add_argument('--electric-library',type=Path,nargs='*',default=[]);p.add_argument('--openems-python');p.add_argument('--openems-dll-directory')
    p.add_argument('--workers',type=int,choices=[1,2],default=1);p.add_argument('--skip-josim',action='store_true')
    p.add_argument('--pitch-um',type=float,default=2.);p.add_argument('--local-pitch-um',type=float,default=.125)
    p.add_argument('--patch-side-um',type=float,default=1.);p.add_argument('--nh',type=int,default=2)
    p.add_argument('--contact-control-pitch-um',type=float,help='Additional coarser local grid after retaining earlier rejected levels')
    p.add_argument('--timeout',type=float,default=900.);p.add_argument('--memory-gib',type=float,default=8.)
    p.add_argument('--dense-threads',type=int,choices=[1,2,3,4],default=2)
    p.add_argument('--executable',default='josim-cli');p.add_argument('--wsl');p.add_argument('--corrector',type=int,default=0)
    p.add_argument('--pump-frequency-hz',type=float,default=5e9);p.add_argument('--pump-power-w',type=float,default=1e-10)
    p.add_argument('--signal-frequency-hz',type=float);p.add_argument('--flux-bias',type=float,default=.05)
    p.add_argument('--flux-pump',type=float,default=0.);p.add_argument('--mixing-order',type=int,choices=[1,2])
    p.add_argument('--output-format',choices=['csv','npy'],default='npy')
    p.add_argument('--electric-port-ohm',type=float,default=1.);p.add_argument('--electric-port-control-ohm',type=float)
    p.add_argument('--electric-timesteps',type=int,default=100000)
    p.add_argument('--josim-step-s',type=float,default=5e-15)
    p.add_argument('--acquisition-start-s',type=float,default=80e-9);p.add_argument('--acquisition-stop-s',type=float,default=160e-9)
    p.add_argument('--port-impedance-ohm',type=float,default=50.)
    p.add_argument('--junction-resistance-ohm',type=float,default=1e4);p.add_argument('--dielectric-resistance-ohm',type=float,default=1e4)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);gds=a.output/'layout.gds'
    if a.gds:
        shutil.copy2(a.gds,gds);shutil.copy2(a.gds.with_suffix('.manifest.json'),gds.with_suffix('.manifest.json'))
    else:generate(JosephsonLayout(family=a.family,cells=a.cells,large_junctions=a.large_junctions),load_process('jj_overlap_illustrative'),gds)
    if json.loads(gds.with_suffix('.manifest.json').read_text())['family']!=a.family:raise ValueError('--family differs from the supplied layout manifest')
    m=read_magnetic(a.magnetic) if a.magnetic else magnetic(gds,a.output,FastHenry.from_installation(a.installation,memory_GiB=a.memory_gib,timeout=a.timeout,threads=a.dense_threads),a.magnetic_primary,a.workers,
            pitch_um=a.pitch_um,local_pitch_um=a.local_pitch_um,patch_side_um=a.patch_side_um,nh=a.nh,contact_pitch_um=a.contact_control_pitch_um)
    em=OpenEMS(a.openems_python,a.openems_dll_directory) if a.openems_python else None
    e=electric(gds,a.output,em,a.electric_library,port_resistance=a.electric_port_ohm,port_control_resistance=a.electric_port_control_ohm,timesteps=a.electric_timesteps)
    circuit,record=back_annotate_layout(gds,e,m,port_impedance=a.port_impedance_ohm,
        junction_resistance=a.junction_resistance_ohm,dielectric_resistance=a.dielectric_resistance_ohm,output=a.output/'backannotation.json')
    if not a.skip_josim:
        from twpa.backends import JoSIM
        compare(circuit,a.output,JoSIM(a.executable,a.wsl,timeout=240,corrector=a.corrector,output_format=a.output_format),
                order=a.mixing_order or (1 if a.flux_pump or a.family in ('rfsquid3wm','snail3wm') else 2),
                pump_frequency=a.pump_frequency_hz,pump_power=a.pump_power_w,signal_frequency=a.signal_frequency_hz,
                flux_bias=a.flux_bias,flux_pump=a.flux_pump,step=a.josim_step_s,start=a.acquisition_start_s,stop=a.acquisition_stop_s)


if __name__=='__main__':main()
