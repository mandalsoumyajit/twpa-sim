"""Generate JPA -> openEMS extraction -> back-annotation -> JoSIM validation.

Run --help for the separate openEMS Python and JoSIM executable paths. Runtime
outputs are local regenerable artifacts. --extract-only omits JoSIM. To reuse
a certified extraction, supply --gds and --extraction together.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
from scipy.constants import h

from twpa import jpa
from twpa.backends import JoSIM
from twpa.jpa_circuit import Tone,Transient,jpa_netlist,ParallelJPA
from twpa.josim_jpa import run_circuit,sideband_measurement
from twpa.layout import load_process
from twpa.layout.generate import generate,JosephsonLayout
from twpa.em import OpenEMS,load_stack
from twpa.em.geometry import capacitor_job
from twpa.em.backannotation import certify,back_annotate_jpa


def db(x): return 10*np.log10(x)


def compare(runner,circuit,root,*,step=2.5e-15,stop=200e-9,start=100e-9,retain_waveforms=False):
    """Compare independent transient RCSJ and reduced single-mode amplifier.

    Uses fixed circuit-derived parameters; no fit to pumped JoSIM gain.
    Reads waveforms before optionally removing only the exact completed CSV.
    """
    transient=Transient(step,stop,start,2e-13)
    r=circuit.reduced_model(); gamma=r.kappa/(2*np.pi)
    def simulate(label,tones):
        result=run_circuit(runner,circuit,root/label,tones,transient)
        if not retain_waveforms:
            # Inputs, logs and run metadata survive. Arrays have been loaded by runner.
            (result.directory/'output.csv').unlink()
        return result
    fp=r.frequency-2*gamma
    fs=r.frequency+gamma*np.array([-.5,0,.5])
    base=simulate('passive-base',[Tone(fp,0)])
    probe=1e-24
    run=simulate('passive-probe',[Tone(f,probe) for f in fs])
    measured=sideband_measurement(run,base,circuit.port_impedance,fp,fs)
    passive=jpa.kerr_response(r,jpa.select_pump_point(jpa.pump_points(r,fp,0)),fs)
    observed=np.asarray(measured['signal_real'])+1j*np.asarray(measured['signal_imag'])
    passive_error=float(max(abs(observed-passive.scattering[:,0,0])))
    n=.28*r.kappa/abs(r.kerr); detuning=2*r.kerr*n
    fp=r.frequency+detuning/(2*np.pi)
    power=n*((r.kappa/2)**2+(detuning-r.kerr*n)**2)*h*fp/r.kappa_external
    point=jpa.select_pump_point(jpa.pump_points(r,fp,power,phase=np.pi/2))
    fs=fp+gamma*np.array([.25,.5])
    tones=[Tone(fp,power)]
    base=simulate('pump-base',tones)
    run=simulate('pump-probe',tones+[Tone(f,probe) for f in fs])
    measured_pump=sideband_measurement(run,base,circuit.port_impedance,fp,fs)
    predicted=jpa.kerr_response(r,point,fs)
    gain_error=db(np.asarray(measured_pump['gain']))-predicted.gain_db
    conversion=abs(predicted.scattering[:,1,0])**2
    idler_error=db(np.asarray(measured_pump['idler_conversion']))-db(conversion)
    passed=bool(passive_error<.02 and max(abs(gain_error))<.25 and max(abs(idler_error))<.25)
    return dict(passed=passed,simulator=runner.probe(),circuit=asdict(circuit),model=asdict(r),
                transient=asdict(transient),probe_power_W=probe,
                passive=dict(measured=measured,max_complex_s11_error=passive_error,
                             model_real=passive.scattering[:,0,0].real.tolist(),model_imag=passive.scattering[:,0,0].imag.tolist()),
                pumped=dict(pump_frequency_Hz=fp,pump_power_W=power,pump_photons=point.photons,
                            measured=measured_pump,model_gain=predicted.gain.tolist(),
                            model_idler_conversion=conversion.tolist(),gain_error_db=gain_error.tolist(),
                            idler_error_db=idler_error.tolist()),
                acceptance=dict(max_complex_s11_error=.02,max_gain_error_db=.25,max_idler_error_db=.25))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--openems-python')
    p.add_argument('--openems-dll-directory')
    p.add_argument('--stack',type=Path)
    p.add_argument('--executable',default='josim-cli')
    p.add_argument('--wsl',help='Optional WSL distribution for JoSIM')
    p.add_argument('--corrector',type=int,default=0)
    p.add_argument('--revision')
    p.add_argument('--gds',type=Path)
    p.add_argument('--extraction',type=Path)
    p.add_argument('--extract-only',action='store_true')
    p.add_argument('--step',type=float,default=2.5e-15)
    p.add_argument('--stop',type=float,default=200e-9)
    p.add_argument('--start',type=float,default=100e-9)
    p.add_argument('--port-impedance',type=float,default=5000.)
    p.add_argument('--retain-waveforms',action='store_true')
    args=p.parse_args()
    if bool(args.gds)!=bool(args.extraction): p.error('--gds and --extraction are required together')
    if not args.gds and not args.openems_python: p.error('--openems-python is required for a new extraction')
    root=args.output.resolve(); root.mkdir(parents=True,exist_ok=False)
    if args.gds:
        gds=args.gds.resolve(); extraction=args.extraction.resolve()
    else:
        gds=root/'jpa.gds'
        generate(JosephsonLayout(),load_process('jj_overlap_illustrative'),gds)
        stack=load_stack(args.stack)
        runner=OpenEMS(args.openems_python,args.openems_dll_directory,timeout=600,threads=1)
        # The accepted primary uses all most-refined controls. Each comparison
        # changes one group only, including both lateral and vertical mesh.
        controls=dict(primary=dict(mesh_um=.5,z_cells=3,margin_um=7,timesteps=100000),
                      mesh=dict(mesh_um=1.,z_cells=2,margin_um=7,timesteps=100000),
                      domain=dict(mesh_um=.5,z_cells=3,margin_um=5,timesteps=100000),
                      duration=dict(mesh_um=.5,z_cells=3,margin_um=7,timesteps=50000))
        results={}
        for name,settings in controls.items():
            print('openEMS',name,settings,flush=True)
            results[name]=runner.run(capacitor_job(gds,stack,**settings),root/name)
            print('  C =',results[name]['capacitance_F'],'F',flush=True)
        certified=certify(results['primary'],[results[k] for k in ('mesh','domain','duration')])
        extraction=root/'extraction-certified.json'
        extraction.write_text(json.dumps(certified,indent=2)+'\n',encoding='utf-8')
        if not certified['convergence']['passed']: raise SystemExit('EM refinement gates failed')
    circuit,annotation=back_annotate_jpa(gds,extraction,port_impedance=args.port_impedance,
                                        output=root/'backannotation.json')
    # Save a reviewable physical circuit even for an extraction-only run.
    (root/'annotated.cir').write_text(jpa_netlist(circuit,[Tone(circuit.reduced_model().frequency,0)],
                                            Transient(args.step,args.stop,args.start,2e-13)),encoding='utf-8')
    print('Back-annotated resonance',circuit.reduced_model().frequency,'Hz',flush=True)
    if args.extract_only: return
    runner=JoSIM(args.executable,args.wsl,timeout=240,corrector=args.corrector,source_revision=args.revision)
    summary=compare(runner,circuit,root,step=args.step,stop=args.stop,start=args.start,
                    retain_waveforms=args.retain_waveforms)
    if args.wsl:
        resolved=subprocess.run(['wsl','-d',args.wsl,'--exec','which',args.executable],
                                capture_output=True,text=True,check=True).stdout.strip()
        binary_hash=subprocess.run(['wsl','-d',args.wsl,'--exec','sha256sum',resolved],
                                   capture_output=True,text=True,check=True).stdout.split()[0]
    else:
        resolved=shutil.which(args.executable) or args.executable
        binary_hash=hashlib.sha256(Path(resolved).read_bytes()).hexdigest()
    summary['simulator_binary_sha256']=binary_hash
    nominal=ParallelJPA(circuit.critical_current,annotation['nominal_shunt_capacitance_F'],
                        circuit.junction_capacitance,circuit.junction_resistance,
                        circuit.port_impedance,circuit.internal_resistance)
    summary['backannotation']=annotation
    summary['nominal_resonance_Hz']=nominal.reduced_model().frequency
    summary['extracted_resonance_Hz']=circuit.reduced_model().frequency
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=summary['passed'],passive_error=summary['passive']['max_complex_s11_error'],
                         gain_error_db=summary['pumped']['gain_error_db'],idler_error_db=summary['pumped']['idler_error_db'])),flush=True)
    if not summary['passed']: raise SystemExit('Back-annotated circuit/model validation failed')


if __name__=='__main__': main()
