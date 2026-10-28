"""Prepare a capacitor coupon from actual exported JPA conductor polygons."""
import hashlib
from pathlib import Path

from ..layout.generate import database, digest
from ..layout.process import _json, Process, positive
from ..layout.validate import validate


def capacitor_job(gds, stack, *, mesh_um=1., margin_um=5., z_cells=2,
                  frequencies_hz=(10e9,20e9,40e9,60e9,80e9,100e9),
                  timesteps=100000, port_resistance=1., crop_margin_um=1.):
    """Extract the shunt-capacitor electrode footprint from a checked JPA GDS.

    The crop is explicit: disconnected leads outside the capacitor coupon are
    not modeled. Film thickness is zero (PEC sheets). Junction capacitance
    is excluded and must be added exactly once by circuit back-annotation.
    """
    import numpy as np
    positive(mesh_um,'mesh_um'); positive(margin_um,'margin_um')
    positive(crop_margin_um,'crop_margin_um'); positive(port_resistance,'port_resistance')
    if isinstance(z_cells,bool) or not isinstance(z_cells,int) or z_cells<2:
        raise ValueError('At least two cells through the capacitor dielectric are required')
    if not isinstance(timesteps,int) or isinstance(timesteps,bool) or timesteps<10000:
        raise ValueError('timesteps must be an integer >=10000')
    fs=np.asarray(frequencies_hz,float)
    if fs.ndim!=1 or len(fs)<3 or np.any(~np.isfinite(fs)) or np.any(fs<=0) or np.any(np.diff(fs)<=0):
        raise ValueError('At least three positive increasing extraction frequencies required')
    gds=Path(gds); manifest_path=gds.with_suffix('.manifest.json')
    m=_json(manifest_path.read_text(encoding='utf-8-sig'))
    if m['family']!='jpa':
        raise ValueError('Initial EM extraction supports single-junction jpa layouts only')
    checked=validate(gds)
    if not checked['passed'] or not checked['gds_matches_manifest']:
        raise ValueError('GDS must pass geometry/DRC and match its manifest before extraction')
    p=Process(m['process']); db=database()
    ly=db.Layout(); ly.read(str(gds)); top=ly.top_cell()
    caps=[c for c in m['components'] if c['kind']=='CAP']
    jjs=[c for c in m['components'] if c['kind']=='JJ']
    if len(caps)!=1 or len(jjs)!=1:
        raise ValueError('Exactly one capacitor and one JJ required')
    cap=caps[0]; bb=db.Box(*(p.coordinate(v) for v in cap['box_um']))
    crop=bb.enlarged(p.coordinate(crop_margin_um))
    center=bb.center(); polygons={}; areas={}
    for role in ('M1','M2'):
        reg=(db.Region(top.begin_shapes_rec(ly.layer(*p.layer(role)))) & db.Region(crop)).merged()
        if reg.count()!=1:
            raise ValueError('Coupon must contain one connected electrode per metal')
        pieces=[]
        for poly in reg.each():
            # Fracture into hole-free trapezoids for CSXCAD; retain actual polygons.
            for trapezoid in poly.decompose_trapezoids():
                pieces.append([[(v.x-center.x)*ly.dbu,(v.y-center.y)*ly.dbu] for v in trapezoid.each_point()])
        polygons[role]=pieces; areas[role]=reg.area()*ly.dbu**2
    marker_um=[(bb.left-center.x)*ly.dbu,(bb.bottom-center.y)*ly.dbu,
               (bb.right-center.x)*ly.dbu,(bb.top-center.y)*ly.dbu]
    width=min(bb.width(),bb.height())*ly.dbu
    port_width=min(2.,width/3)
    if port_width<mesh_um:
        raise ValueError('Mesh too coarse to resolve the port and capacitor interior')
    junction_area=db.Region(top.begin_shapes_rec(ly.layer(*p.layer('JJ')))).area()*ly.dbu**2
    junction=dict(drawn_area_um2=junction_area,fabricated_area_assumed_um2=junction_area,
                  Ic_A=junction_area*p.parameter('critical_current_density','uA/um2')*1e-6,
                  Cj_F=junction_area*p.parameter('junction_specific_capacitance','fF/um2')*1e-15)
    source=dict(gds_sha256=hashlib.sha256(gds.read_bytes()).hexdigest(),
                process_sha256=m['process_sha256'],design_sha256=m['design_sha256'],
                manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                process_id=p.id,component_id=cap['id'],junction_id=jjs[0]['id'],
                junction_estimates=junction,nominal_capacitance_F=cap['estimates']['C_F'])
    return dict(schema_version=1,kind='jpa_mim_capacitor',source=source,
                stack=stack.to_dict(),stack_sha256=digest(stack.to_dict()),
                geometry=dict(polygons_um=polygons,electrode_areas_um2=areas,
                              marker_box_um=marker_um,crop_margin_um=crop_margin_um,
                              port_width_um=port_width),
                settings=dict(mesh_um=mesh_um,margin_um=margin_um,z_cells=z_cells,
                              frequencies_hz=fs.tolist(),timesteps=timesteps,port_resistance=port_resistance,
                              excitation_cutoff_hz=2e12),
                scope='Component-level PEC MIM coupon from GDS. Interconnect, package, bias structures and junction barrier excluded.')
