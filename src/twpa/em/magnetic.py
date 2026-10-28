"""GDS film-current mesh and fixed finite electrodes for London FastHenry.

Method follows the finite-patch/explicit-return workflow reviewed in the local
spec2si-seeqc project. This independent adapter imports no vendor process data.
"""
from dataclasses import asdict,dataclass
from importlib.resources import files
import math
from pathlib import Path
import numpy as np

from ..layout.generate import database,digest
from ..layout.process import _json,positive,integer
from ..layout.validate import validate


@dataclass(frozen=True)
class MagneticStack:
    id: str = 'illustrative_london_v1'
    status: str = 'illustrative'
    temperature_K: float = .02
    m1_thickness_um: float = .2
    m2_thickness_um: float = .2
    dielectric_gap_um: float = .0664
    m1_lambda_um: float = .09
    m2_lambda_um: float = .09
    via_lambda_um: float = .09

    def __post_init__(self):
        if self.status not in ('illustrative','local_unqualified') or not self.id:
            raise ValueError('Explicitly unqualified magnetic stack required')
        for name,value in asdict(self).items():
            if name not in ('id','status'): positive(value,name)

    def film(self,role):
        if role=='M1': return (-self.m1_thickness_um/2,self.m1_thickness_um,self.m1_lambda_um)
        if role=='M2': return (self.dielectric_gap_um+self.m2_thickness_um/2,self.m2_thickness_um,self.m2_lambda_um)
        raise ValueError('Only M1/M2 supported')


def load_magnetic_stack(path=None):
    data=_json(Path(path).read_text(encoding='utf-8') if path else files('twpa.em').joinpath('stacks','london_illustrative.json').read_text(encoding='utf-8'))
    if data.pop('schema_version',None)!=1:raise ValueError('Unsupported magnetic stack')
    try:return MagneticStack(**data)
    except TypeError as e:raise ValueError(f'Invalid magnetic stack fields: {e}') from e


def _cuts(values,pitch,bands):
    anchors=sorted(set(round(float(v),9) for v in values))
    out=[anchors[0]]
    for a,b in zip(anchors,anchors[1:]):
        step=min([pitch]+[q for lo,hi,q in bands if lo-1e-9<=(a+b)/2<=hi+1e-9])
        n=max(1,math.ceil((b-a)/step-1e-9))
        out.extend(round(a+(b-a)*i/n,9) for i in range(1,n+1))
    return out


def jpa_magnetic_job(gds,stack=None,*,pitch_um=3.,local_pitch_um=.5,nh=3,
                     patch_side_um=1.,lambda_scale=1.,frequencies_hz=(1e6,1e7,1e8,1e9,1e10,1e11)):
    """Mesh ALL M1/M2 conductors, filled contacts and an optional pickup loop.

    Electrode footprints are held fixed when refining the current basis.
    Physical ports run within each connected conductor; a balanced transform
    forms the closed tank and pickup-current loops. No ground plane is idealized.
    """
    import hashlib
    stack=stack or MagneticStack()
    for n,v in [('pitch_um',pitch_um),('local_pitch_um',local_pitch_um),('patch_side_um',patch_side_um)]: positive(v,n)
    integer(nh,'nh')
    if isinstance(lambda_scale,bool) or not math.isfinite(lambda_scale) or lambda_scale<=0: raise ValueError('Positive lambda_scale required; use an explicitly converged small-lambda PEC limit')
    fs=np.asarray(frequencies_hz,float)
    if fs.ndim!=1 or len(fs)<1 or np.any(~np.isfinite(fs)) or np.any(fs<=0) or np.any(np.diff(fs)<=0):
        raise ValueError('Increasing positive finite frequencies required')
    path=Path(gds); m=_json(path.with_suffix('.manifest.json').read_text(encoding='utf-8'))
    if m['family']!='jpa': raise ValueError('Magnetic extraction currently supports single-JJ JPAs')
    checked=validate(path)
    if not checked['passed'] or not checked['gds_matches_manifest']: raise ValueError('Layout must pass DRC and match its manifest')
    db=database(); ly=db.Layout(); ly.read(str(path)); top=ly.top_cell()
    from ..layout import Process
    process=Process(m['process']); regions={role:db.Region(top.begin_shapes_rec(ly.layer(*process.layer(role)))).merged() for role in ('M1','M2','VIA')}
    rects={};xs=[];ys=[]
    for role in ('M1','M2'):
        rects[role]=[]
        for poly in regions[role].each():
            for shape in poly.decompose_trapezoids():
                if not shape.is_box(): raise ValueError('Current magnetic adapter requires Manhattan geometry')
                b=shape.bbox(); r=[v*ly.dbu for v in (b.left,b.bottom,b.right,b.top)]
                rects[role].append(r);xs+=r[::2];ys+=r[1::2]
    components={c['kind']:c for c in m['components'] if c['kind'] in ('JJ','CAP')}
    if set(components)!= {'JJ','CAP'}: raise ValueError('Single JJ and shunt capacitor required')
    patches=[]
    for kind in ('JJ','CAP'):
        b=components[kind]['box_um'];x=(b[0]+b[2])/2;y=(b[1]+b[3])/2
        if patch_side_um>min(b[2]-b[0],b[3]-b[1]):raise ValueError('Electrode patch exceeds device window')
        for role in ('M1','M2'): patches.append(dict(name=kind+'_'+role,role=role,box=[x-patch_side_um/2,y-patch_side_um/2,x+patch_side_um/2,y+patch_side_um/2]))
    pins={p['name']:p for p in m['pins']}
    coupled=m['design'].get('inductive_coupler',False)
    if coupled:
        for name in ('coil_lo','coil_hi'):
            x,y=pins[name]['at_um'];patches.append(dict(name=name,role='M1',box=[x-patch_side_um/2,y-patch_side_um/2,x+patch_side_um/2,y+patch_side_um/2]))
    else:
        for name,pin in [('EXT','in'),('GROUND','ground_rail')]:
            x,y=pins[pin]['at_um'];patches.append(dict(name=name,role='M1',box=[x-patch_side_um/2,y-patch_side_um/2,x+patch_side_um/2,y+patch_side_um/2]))
    vias=[]
    for poly in regions['VIA'].each():
        if not poly.is_box():raise ValueError('Rectangular contacts required')
        b=poly.bbox();vias.append([v*ly.dbu for v in (b.left,b.bottom,b.right,b.top)])
    for patch in patches:
        box=db.Box(*(process.coordinate(v) for v in patch['box']))
        if not (db.Region(box)-regions[patch['role']]).is_empty():raise ValueError('Finite electrode not fully covered by actual metal')
    xbands=[];ybands=[]
    for b in [p['box'] for p in patches]+vias:
        xs+=b[::2];ys+=b[1::2];xbands.append((b[0],b[2],local_pitch_um));ybands.append((b[1],b[3],local_pitch_um))
    xe=_cuts(xs,pitch_um,xbands);ye=_cuts(ys,pitch_um,ybands)
    if (len(xe)-1)*(len(ye)-1)>500000:raise ValueError('Current-grid candidate count exceeds the supported local extraction envelope')
    nodes={};coords={};segments=[];equivs=[];edge_pairs={};lines=['* twpa-sim full JPA current basis','.units um','.default sigma=0']
    def node(role,x,y,z):
        key=(role,round(x,9),round(y,9),round(z,9))
        if key not in nodes:
            n='N'+str(len(nodes));nodes[key]=n;coords[n]=[x,y,z]
            lines.append(n+f' x={x:.12g} y={y:.12g} z={z:.12g}')
        return nodes[key]
    def edge(a,b,width,height,lam,axis,role,vertical=False):
        if np.linalg.norm(np.array(coords[a])-coords[b])<1e-10:raise ValueError('Zero-length current branch')
        vec=(0,1,0) if axis=='x' else (1,0,0)
        pair=tuple(sorted((a,b)));signature=(round(width,9),round(height,9),round(lam,9),axis,vertical)
        if pair in edge_pairs:
            if edge_pairs[pair]!=signature:raise ValueError('Inconsistent overlapping physical current branch')
            return
        edge_pairs[pair]=signature
        k='E'+str(len(segments));nheight=1 if vertical else nh
        lines.append(f'{k} {a} {b} w={width:.12g} h={height:.12g} lambda={lam*lambda_scale:.12g} sigma=0 nwinc=1 nhinc={nheight} rw=1 rh=1 wx={vec[0]} wy={vec[1]} wz={vec[2]}')
        segments.append(dict(a=a,b=b,w_um=width,h_um=height,lambda_um=lam*lambda_scale,axis=axis,role=role,nhinc=nheight))
    tiles={};areas={}
    for role in ('M1','M2'):
        z,t,lam=stack.film(role);cells={}
        for i,(a,b) in enumerate(zip(xe,xe[1:])):
            x=(a+b)/2
            for j,(c,d) in enumerate(zip(ye,ye[1:])):
                y=(c+d)/2
                if any(r[0]<x<r[2] and r[1]<y<r[3] for r in rects[role]):
                    cells[i,j]=([a,c,b,d],node(role,x,y,z))
        areas[role]=sum((b[2]-b[0])*(b[3]-b[1]) for b,n in cells.values())
        if not math.isclose(areas[role],regions[role].area()*ly.dbu**2,abs_tol=1e-6):raise ValueError('Mesh does not cover exact conductor area')
        for (i,j),(b,center) in cells.items():
            x,y=(b[0]+b[2])/2,(b[1]+b[3])/2
            for ij,point,axis,w in [((i-1,j),(b[0],y),'x',b[3]-b[1]),((i+1,j),(b[2],y),'x',b[3]-b[1]),
                                     ((i,j-1),(x,b[1]),'y',b[2]-b[0]),((i,j+1),(x,b[3]),'y',b[2]-b[0])]:
                if ij in cells:edge(center,node(role,*point,z),w,t,lam,axis,role)
        tiles[role]=cells
    terminals={}
    for patch in patches:
        role=patch['role'];z,t,lam=stack.film(role);b=patch['box'];selected=[];area=0
        # Facing surfaces: M1 top and M2 bottom. Fixed patch feed itself is excluded.
        surface=z+t/2 if role=='M1' else z-t/2
        for r,center in tiles[role].values():
            x,y=(r[0]+r[2])/2,(r[1]+r[3])/2
            if b[0]<x<b[2] and b[1]<y<b[3]:
                endpoint=node(role,x,y,surface);edge(center,endpoint,r[2]-r[0],r[3]-r[1],lam,'z','terminal_half_film',True)
                selected.append(endpoint);area+=(r[2]-r[0])*(r[3]-r[1])
        if not selected or not math.isclose(area,patch_side_um**2,abs_tol=1e-7):raise ValueError('Incomplete electrode footprint')
        terminals[patch['name']]=selected[0]
        for n in selected[1:]:lines.append('.equiv '+selected[0]+' '+n);equivs.append((selected[0],n))
    for b in vias:
        area=0
        for ij,(r,lower) in tiles['M1'].items():
            x,y=(r[0]+r[2])/2,(r[1]+r[3])/2
            if b[0]<x<b[2] and b[1]<y<b[3]:
                upper=tiles['M2'][ij][1];z1,t1,l1=stack.film('M1');z2,t2,l2=stack.film('M2')
                chain=[lower,node('M1',x,y,z1+t1/2),node('M2',x,y,z2-t2/2),upper]
                for a,c,lam in zip(chain,chain[1:],[l1,stack.via_lambda_um,l2]):edge(a,c,r[2]-r[0],r[3]-r[1],lam,'z','via_film_or_fill',True)
                area+=(r[2]-r[0])*(r[3]-r[1])
        if not math.isclose(area,(b[2]-b[0])*(b[3]-b[1]),abs_tol=1e-7):raise ValueError('Incomplete via fill')
    ports=[dict(name='signal_path',positive=terminals['CAP_M2'],negative=terminals['JJ_M1'],positive_patch='CAP_M2',negative_patch='JJ_M1'),
           dict(name='return_path',positive=terminals['JJ_M2'],negative=terminals['CAP_M1'],positive_patch='JJ_M2',negative_patch='CAP_M1')]
    if coupled:ports.append(dict(name='pickup',positive=terminals['coil_lo'],negative=terminals['coil_hi'],positive_patch='coil_lo',negative_patch='coil_hi'))
    else:
        definitions=[('jj_signal','EXT','JJ_M1'),('cap_signal','EXT','CAP_M2'),
                     ('jj_return','JJ_M2','GROUND'),('cap_return','CAP_M1','GROUND')]
        ports=[dict(name=name,positive=terminals[a],negative=terminals[b],positive_patch=a,negative_patch=b) for name,a,b in definitions]
    # Audit every physical source/sink pair in the full conductive graph.
    parent={n:n for n in coords}
    def find(n):
        while parent[n]!=n:parent[n]=parent[parent[n]];n=parent[n]
        return n
    for a,b in [(e['a'],e['b']) for e in segments]+equivs:parent[find(a)]=find(b)
    for port in ports:
        if find(port['positive'])!=find(port['negative']):raise ValueError('Disconnected FastHenry physical port')
        lines.append('.external '+port['positive']+' '+port['negative'])
    if sum(e['nhinc'] for e in segments)>20000:raise ValueError('Current basis exceeds the declared 20,000-filament local envelope')
    decades=math.log10(fs[-1]/fs[0]);ndec=1 if len(fs)==1 else (len(fs)-1)/decades
    expected=np.geomspace(fs[0],fs[-1],round(decades*ndec)+1)
    if len(expected)!=len(fs) or not np.allclose(expected,fs,rtol=1e-10):raise ValueError('FastHenry frequency grid must be logarithmic')
    lines += [f'.freq fmin={fs[0]:.12g} fmax={fs[-1]:.12g} ndec={ndec:.12g}','.end','']
    transform=[[1.,0.],[1.,0.],[0.,1.]] if coupled else [[1.,0.],[0.,1.],[1.,0.],[0.,1.]]
    return dict(schema_version=1,kind='jpa_london_network',deck='\n'.join(lines),
                source={k:m[k] for k in ('gds_sha256','process_sha256','design_sha256')},
                manifest_sha256=hashlib.sha256(path.with_suffix('.manifest.json').read_bytes()).hexdigest(),
                stack=asdict(stack),settings=dict(pitch_um=pitch_um,local_pitch_um=local_pitch_um,nh=nh,
                  patch_side_um=patch_side_um,lambda_scale=lambda_scale,frequencies_hz=fs.tolist()),
                ports=ports,reduction=transform,patches=patches,topology='pickup' if coupled else 'galvanic',
                audit=dict(conductor_areas_um2=areas,conductor_components=len({find(n) for n in coords}),
                           nodes=len(coords),segments=len(segments),filaments=sum(e['nhinc'] for e in segments),
                           vias_um=vias,physical_ports_connected=True,ideal_equivalence_scope='within each declared fixed local electrode only'),
                segments=segments,coordinates=coords,
                scope='Full actual M1/M2 finite-film current grid and filled contacts; London MQS. JJ and shunt electric energy are separate circuit elements. Local electrode feeds excluded.')
