"""Standalone openEMS worker. No twpa or KLayout imports in this interpreter.

This file is copied into each immutable run directory for reproduction. The
caller supplies JSON geometry extracted by KLayout and an explicit stack.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time


def runtime(dll_directory=None):
    handles=[]
    if dll_directory and os.name=='nt':
        handles.append(os.add_dll_directory(str(Path(dll_directory).resolve())))
    import openEMS
    import CSXCAD
    return openEMS,CSXCAD,handles


def run(job,root,threads,dll_directory=None):
    import numpy as np
    em,csx,handles=runtime(dll_directory)
    stack=job['stack']; settings=job['settings']; geo=job['geometry']
    dx=settings['mesh_um']; h=stack['dielectric_thickness_um']
    margin=settings['margin_um']; zc=settings['z_cells']
    fdtd=em.openEMS(NrTS=settings['timesteps'],EndCriteria=1e-5)
    fdtd.SetGaussExcite(0,settings['excitation_cutoff_hz'])
    fdtd.SetBoundaryCond(['MUR']*6)
    c=csx.ContinuousStructure(); fdtd.SetCSX(c)
    grid=c.GetGrid(); grid.SetDeltaUnit(1e-6)
    all_points=np.array([xy for shapes in geo['polygons_um'].values() for poly in shapes for xy in poly])
    low=np.min(all_points,axis=0)-margin; high=np.max(all_points,axis=0)+margin
    width=geo['port_width_um']; half=width/2
    mesh=[]
    for axis in range(2):
        vertices=np.unique(all_points[:,axis])
        # All polygon vertices and port edges are exact grid planes.
        lines=np.unique(np.r_[np.arange(low[axis],high[axis]+dx/2,dx),low[axis],high[axis],vertices,-half,half,0])
        if np.min(np.diff(lines))<dx*1e-5:
            raise ValueError('Near-duplicate mesh lines; choose a compatible mesh step')
        grid.AddLine('xy'[axis],lines); mesh.append(lines)
    sub=stack['substrate_thickness_um']
    zs=np.unique(np.r_[-sub,-sub*.66,-sub*.33,-.5,-.2,-h,-h/2,
                       np.linspace(0,h,zc+1),h*1.5,h*2,h+.2,h+.5,h+1,h+2,h+sub])
    # Restrict substrate lines to the declared volume; keep air above M2.
    zs=zs[zs>=-sub]
    if np.min(np.diff(zs))<1e-8:
        raise ValueError('Near-duplicate vertical mesh planes')
    grid.AddLine('z',zs); mesh.append(zs)
    dielectric=c.AddMaterial('dielectric',epsilon=stack['dielectric_epsilon_r'])
    dielectric.AddBox([*low,0],[*high,h])
    substrate=c.AddMaterial('substrate',epsilon=stack['substrate_epsilon_r'])
    substrate.AddBox([*low,-sub],[*high,0])
    for role,z in (('M1',0),('M2',h)):
        metal=c.AddMetal(role)
        for polygon in geo['polygons_um'][role]:
            metal.AddPolygon(np.asarray(polygon).T.tolist(),norm_dir='z',elevation=z,priority=10)
    port=fdtd.AddLumpedPort(1,settings['port_resistance'],[-half,-half,0],[half,half,h],
                           'z',1,priority=20)
    c.Write2XML(str(root/'model.xml'))
    started=time.monotonic()
    fdtd.Run(str(root),verbose=0,cleanup=False,numThreads=threads)
    freq=np.asarray(settings['frequencies_hz'])
    port.CalcPort(str(root),freq)
    admittance=port.if_tot/port.uf_tot  # openEMS uses exp(+i omega t).
    cs=admittance.imag/(2*np.pi*freq)
    fitted=float(np.mean(cs))
    residual=float(np.max(abs(admittance-1j*2*np.pi*freq*fitted)/abs(admittance)))
    u=np.asarray(port.ut_tot); current=np.asarray(port.it_tot)
    tail=max(float(np.max(abs(u[int(.9*len(u)):]))/np.max(abs(u))),
             float(np.max(abs(current[int(.9*len(current)):]))/np.max(abs(current))))
    # Passive conductance is a numerical/finite-domain diagnostic, not a fitted
    # internal resistor or an inferred superconducting/material loss.
    negative=float(np.max(np.maximum(-admittance.real,0)/abs(admittance)))
    accepted=bool(np.isfinite(fitted) and fitted>0 and residual<.02 and tail<1e-3 and negative<.01)
    result=dict(schema_version=1,kind=job['kind'],source=job['source'],stack=stack,
                stack_sha256=job['stack_sha256'],geometry=geo,settings=settings,
                simulator=dict(name='openEMS',version=em.__version__,csxcad_version=csx.__version__),
                mesh=dict(lines=[len(x) for x in mesh],cells=int(np.prod([len(x)-1 for x in mesh])),
                          min_steps_um=[float(np.min(np.diff(x))) for x in mesh]),
                elapsed_seconds=time.monotonic()-started,frequency_hz=freq.tolist(),
                admittance_real_S=admittance.real.tolist(),admittance_imag_S=admittance.imag.tolist(),
                capacitance_by_frequency_F=cs.tolist(),capacitance_F=fitted,
                parallel_C_fit_max_relative_residual=residual,port_tail_relative=tail,
                negative_conductance_relative=negative,fit_accepted=accepted,
                convergence='Unestablished until independent mesh/domain/time refinements are compared.',
                convention='openEMS exp(+i omega t); downstream twpa uses exp(-i omega t). C is real.',
                scope=job['scope'])
    (root/'extraction.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path)
    parser.add_argument('--dll-directory')
    parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--probe',action='store_true')
    args=parser.parse_args()
    if args.probe:
        em,csx,handles=runtime(args.dll_directory)
        print('TWPA_OPENEMS_JSON '+json.dumps(dict(openems=em.__version__,csxcad=csx.__version__)))
        return
    if args.threads<1 or args.input is None:
        raise ValueError('Input job and positive thread count required')
    job=json.loads(args.input.read_text(encoding='utf-8'))
    result=run(job,args.input.resolve().parent,args.threads,args.dll_directory)
    print('TWPA_OPENEMS_JSON '+json.dumps(dict(fit_accepted=result['fit_accepted'],capacitance_F=result['capacitance_F'])))


if __name__=='__main__':
    main()
