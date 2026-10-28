"""Parameterized research-layout generators; geometry is in micrometres."""
from dataclasses import asdict, dataclass, fields
import hashlib
import json
from pathlib import Path

from .process import Process, integer, positive, _json, _keys


def database():
    try:
        import klayout.db as db
        return db
    except ImportError as e:
        raise ImportError('Install twpa-sim[layout] for layout generation/checking') from e


@dataclass(frozen=True)
class KITLayout:
    """Serpentine IDC-loaded KIT. Loading is symmetric within each supercell."""
    rows: int = 16
    periods_per_row: int = 55
    unloaded_cells: int = 48
    loaded_cells: int = 6
    turn_periods: int = 3
    pitch: float = 1.5
    spine_width: float = .580
    spine_gap: float = .500
    finger_width: float = .500
    finger_gap: float = .250
    finger_length: float = 75.500
    loaded_finger_length: float = 27.500
    turn_radius: float = 154.
    bias_ratio: float = .20
    taper_cells: int = 6
    arc_points: int = 1025
    launch_width: float = 103.418
    launch_gap: float = 24.581
    launch_length: float = 155.
    taper_length: float = 88.
    die_size: float = 5000.
    pad_width: float = 90.
    pad_height: float = 65.

    @classmethod
    def from_process(cls, process, **overrides):
        """Seed geometry from reference parameters, then apply explicit overrides."""
        names={f.name for f in fields(cls)}
        seed={name:process.parameter(name,'um') for name in process.to_dict()['parameters']
              if name in names}
        return cls(**(seed | overrides))

    def validate(self, p):
        for f in fields(self):
            v = getattr(self, f.name)
            if f.name in ('rows', 'periods_per_row', 'unloaded_cells', 'loaded_cells', 'turn_periods', 'taper_cells', 'arc_points'):
                integer(v, f.name)
            elif f.name == 'bias_ratio':
                if not isinstance(v, (int, float)) or not 0 <= v <= 1:
                    raise ValueError('bias_ratio must be between 0 and 1')
            else:
                positive(v, f.name); p.coordinate(v)
        if self.rows < 2 or self.unloaded_cells % 2 or self.arc_points < 33:
            raise ValueError('At least two rows, even unloaded_cells and >=33 arc_points required')
        if p.grid != p.dbu:
            raise ValueError('Curved KIT generator currently requires grid_um == dbu_um')
        if self.pitch <= self.finger_width + 2*self.finger_gap:
            raise ValueError('Finger slots leave no inter-cell ground finger')
        if self.loaded_finger_length > self.finger_length or self.turn_radius <= self.finger_length + 5:
            raise ValueError('Loaded fingers or turn radius are inconsistent')
        if self.launch_width/2 + self.launch_gap >= self.turn_radius:
            raise ValueError('Launch exceeds ground shoulder')
        if self.pad_width > self.launch_length - 2 or self.pad_height > self.launch_width - 2:
            raise ValueError('Signal bump pad does not fit launch')
        if self.taper_cells*2 >= (self.unloaded_cells+self.loaded_cells)*self.turn_periods:
            raise ValueError('Turn tapers occupy the entire turn')
        p.require_roles(dict(DIE='boundary',M1='conductor',GND_PAD='conductor',
                             SIG_PAD='conductor',OXIDE='opening'))


@dataclass(frozen=True)
class JosephsonLayout:
    """Illustrative overlap-JJ devices with MIM shunts and explicit circuit intent.

    family: jpa, squid_jpa, jtwpa4wm, rfsquid3wm, or snail3wm.
    A JPA is a one-port lumped LC device; coupling/bias circuitry is external.
    """
    family: str = 'jpa'
    cells: int = 1
    large_junctions: int = 3
    junction_width: float = 2.
    small_junction_width: float = 1.
    capacitor_width: float = 12.
    capacitor_height: float = 12.
    wire_width: float = 2.
    via_width: float = 1.
    cell_pitch: float = 60.
    branch_height: float = 20.
    die_margin: float = 20.

    def validate(self, p):
        if self.family not in ('jpa', 'squid_jpa', 'jtwpa4wm', 'rfsquid3wm', 'snail3wm'):
            raise ValueError('Unknown Josephson layout family')
        integer(self.cells, 'cells'); integer(self.large_junctions, 'large_junctions', 2)
        if self.family in ('jpa', 'squid_jpa') and self.cells != 1:
            raise ValueError('A resonant JPA uses cells=1')
        for f in fields(self):
            if f.name in ('family', 'cells', 'large_junctions'):
                continue
            v = getattr(self, f.name); positive(v, f.name); p.coordinate(v)
        if not .5 <= self.junction_width <= 4:
            raise ValueError('This primitive supports JJ widths .5 to 4 um')
        if self.family=='snail3wm' and not .5 <= self.small_junction_width <= self.junction_width:
            raise ValueError('SNAIL requires .5 <= small width <= large width')
        if not 1 <= self.wire_width <= 3 or not .5 <= self.via_width <= 1.2:
            raise ValueError('Wire/contact dimensions exceed primitive support')
        if self.capacitor_width < 2 or self.capacitor_width > self.cell_pitch/3 or self.capacitor_height < 2:
            raise ValueError('Capacitor dimensions exceed primitive support')
        if self.cell_pitch < 40 or self.branch_height < 16 or self.die_margin < 6:
            raise ValueError('Insufficient routing/branch/die clearance')
        p.require_roles(dict(DIE='boundary',M1='conductor',M2='conductor',VIA='opening',
                             JJ='junction',CAP='capacitor',IND='inductor'))


class Canvas:
    def __init__(self, p):
        self.db = database(); self.p = p; self.ly = self.db.Layout(); self.ly.dbu = p.dbu
        self.top = self.ly.create_cell('TWPA'); self.components = []; self.pins = []; self.notes = []

    def insert(self, role, polygon):
        self.top.shapes(self.ly.layer(*self.p.layer(role))).insert(polygon)

    def box(self, role, x1, y1, x2, y2):
        if x2 <= x1 or y2 <= y1:
            raise ValueError('Degenerate rectangle')
        box = self.db.Box(*(self.p.coordinate(v) for v in (x1,y1,x2,y2)))
        self.insert(role, box)
        return box

    def wire(self, points, width):
        self.insert('M1', self.db.Path([self.db.Point(self.p.coordinate(x),self.p.coordinate(y)) for x,y in points], self.p.coordinate(width)))

    def pin(self, name, net, x, y, role='M1'):
        self.pins.append(dict(name=name, net=net, layer=role, at_um=[x,y]))

    def component(self, kind, a, b, x1, y1, x2, y2, **estimates):
        self.box(kind, x1, y1, x2, y2)
        self.components.append(dict(id=f'{kind}{len(self.components)+1}', kind=kind,
                                    nets=[a,b], box_um=[x1,y1,x2,y2], estimates=estimates))

    def jj(self, x, y, a, b, width, wire_width=2, via_width=1):
        """Left M1 electrode -> tunnel window -> M2 bridge -> contact -> right M1."""
        self.box('M1', x-8, y-3, x-2+width+.5, y+3)
        self.box('M1', x+4, y-3, x+8, y+3)
        self.box('M2', x-2.5, y-width/2-.5, x+6.5, y+width/2+.5)
        self.box('M2', x+4.5, y-1.5, x+7.5, y+1.5)
        self.box('VIA', x+6-via_width/2, y-via_width/2, x+6+via_width/2, y+via_width/2)
        area = width**2
        self.component('JJ', a,b, x-2,y-width/2,x-2+width,y+width/2,
                       drawn_area_um2=area, fabricated_area_assumed_um2=area,
                       Ic_A=area*self.p.parameter('critical_current_density','uA/um2')*1e-6,
                       Cj_F=area*self.p.parameter('junction_specific_capacitance','fF/um2')*1e-15)

    def capacitor(self, x, net, width, height, wire_width, via_width, ground_y):
        # Lower electrode connects to ground; upper electrode contacts the signal node.
        y = -12-height/2
        self.box('M1', x-width/2-1,y-height/2-1,x+width/2+1,y+height/2+1)
        self.wire([(x,y),(x,ground_y)], wire_width)
        self.box('M2', x-width/2-.5,y-height/2-.5,x+width/2+.5,y+height/2+.5)
        self.box('M2', x-1.5,y,x+1.5,1.5)
        self.box('M1', x-1.5,-1.5,x+1.5,1.5)
        self.box('VIA', x-via_width/2,-via_width/2,x+via_width/2,via_width/2)
        area=width*height
        self.component('CAP',net,'gnd',x-width/2,y-height/2,x+width/2,y+height/2,
                       area_um2=area, C_F=area*self.p.parameter('capacitor_specific_capacitance','fF/um2')*1e-15)


def _josephson(c, d):
    w=d.wire_width; pitch=d.cell_pitch
    ground_y=-24-d.capacitor_height
    nlarge = d.large_junctions if d.family == 'snail3wm' else 1
    span=pitch*nlarge
    for i in range(d.cells):
        x=i*span; end=x+span
        a=f'n{i}'; b='gnd' if d.family in ('jpa','squid_jpa') else f'n{i+1}'
        # One junction on the lower branch, with separately defined upper branch.
        mid=x+span/2
        width=d.small_junction_width if d.family=='snail3wm' else d.junction_width
        c.wire([(x,0),(mid-8,0)],w); c.jj(mid,0,a,b,width,w,d.via_width)
        c.wire([(mid+8,0),(end,0)],w)
        if d.family in ('squid_jpa','rfsquid3wm','snail3wm'):
            h=d.branch_height
            c.wire([(x,0),(x,h),(x+pitch/2-8,h)],w)
            if d.family == 'rfsquid3wm':
                # A wire loop branch: IND marks a cross-section for connectivity extraction.
                c.wire([(x+pitch/2-8,h),(end,h),(end,0)],w)
                c.component('IND',a,b,x+pitch/2-2,h-w/2-.5,x+pitch/2+2,h+w/2+.5,
                            inductance_H=None, requires_em_extraction=True)
            else:
                for j in range(nlarge):
                    xx=x+(j+.5)*pitch
                    aa=a if j==0 else f'loop{i}_{j}'
                    bb=b if j==nlarge-1 else f'loop{i}_{j+1}'
                    c.jj(xx,h,aa,bb,d.junction_width,w,d.via_width)
                    if j<nlarge-1:
                        c.wire([(xx+8,h),(xx+pitch-8,h)],w)
                c.wire([(end-pitch/2+8,h),(end,h),(end,0)],w)
        if i==0:
            c.pin('in',a,x,0)
        if d.family in ('jpa','squid_jpa'):
            c.wire([(end,0),(end,ground_y)],w)
        factor=.5 if i==0 and d.family not in ('jpa','squid_jpa') else 1
        c.capacitor(x,a,d.capacitor_width,d.capacitor_height*factor,w,d.via_width,ground_y)
        if i==d.cells-1:
            c.pin('out' if b!='gnd' else 'ground',b,end,0)
    if d.family not in ('jpa','squid_jpa'):
        c.capacitor(d.cells*span,f'n{d.cells}',d.capacitor_width,d.capacitor_height/2,w,d.via_width,ground_y)
    c.wire([(-12,ground_y),(d.cells*span+12,ground_y)],w)
    c.pin('ground_rail','gnd',0,ground_y)
    bbox=c.top.bbox().enlarged(c.p.coordinate(d.die_margin))
    c.insert('DIE',bbox)
    c.notes += ['Illustrative overlap stack; ports are galvanic pads/terminals with external coupling and bias circuitry.',
                'Electrical estimates assume fabricated area equals drawn area; no lithographic bias is applied.',
                'MIM capacitance uses area only. Inductance, fringing, lead capacitance and flux coupling require EM extraction.',
                'Traveling-wave layouts use symmetric half-shunts at both end nodes, matching the finite-chain circuit convention.']


def _kit(c,d):
    from .kit import primitives
    db=c.db
    sc, right, left, launch_left, launch_right, pitch, nc, radius, wm=primitives(c.ly,d)
    length=d.periods_per_row*nc*pitch
    xl=-length/2; xr=length/2
    ys=[(d.rows-1-2*r)*radius for r in range(d.rows)]
    # Preserve the original comb topology and constant-width radial turn fingers.
    for y in ys:
        for j in range(d.periods_per_row):
            t=db.Trans(db.Vector(c.p.coordinate(xl+j*nc*pitch), c.p.coordinate(y)))
            for poly in sc:
                c.insert('M1',poly.transformed(t))
    for i in range(d.rows-1):
        t=db.Trans(db.Vector(c.p.coordinate(xr if i%2==0 else xl),c.p.coordinate((ys[i]+ys[i+1])/2)))
        for poly in right if i%2==0 else left:
            c.insert('M1',poly.transformed(t))
    bay=d.launch_length+d.taper_length
    for i,y in ((0,ys[0]),(d.rows-1,ys[-1])):
        is_left=i%2==1 or i==0
        x=xl if is_left else xr
        t=db.Trans(db.Vector(c.p.coordinate(x),c.p.coordinate(y)))
        for poly in launch_left if is_left else launch_right:
            c.insert('M1',poly.transformed(t))
        xp=x+(-1 if is_left else 1)*(d.taper_length+d.launch_length/2)
        c.box('SIG_PAD',xp-d.pad_width/2,y-d.pad_height/2,xp+d.pad_width/2,y+d.pad_height/2)
        c.pin('in' if i==0 else 'out','signal',xp,y)
    # Ground landing pads between straight rows; all are checked against the actual metal.
    pw=min(d.pad_width,length/3)
    for i in range(d.rows-1):
        y=(ys[i]+ys[i+1])/2
        c.box('GND_PAD',-pw/2,y-d.pad_height/2,pw/2,y+d.pad_height/2)
    # Separate outer/inner ground components are intentional in this single-layer design.
    c.pin('ground_outer','gnd_outer',0,ys[0]+radius-10)
    c.pin('ground_inner','gnd_inner',0,(ys[0]+ys[1])/2)
    bb=c.top.bbox()
    margin=10
    c.box('OXIDE',bb.left*c.p.dbu-margin,bb.bottom*c.p.dbu-margin,
          bb.right*c.p.dbu+margin,bb.top*c.p.dbu+margin)
    half=d.die_size/2
    c.box('DIE',-half,-half,half,half)
    c.notes += [f'{d.rows*d.periods_per_row*nc+(d.rows-1)*d.turn_periods*nc} distributed cells; turn maximum spine width {wm:.9g} um.',
                'Turn matching preserves the original biased sheet-squares approximation; EM impedance is not validated.',
                'Two disconnected ground combs require package grounding. Serpentine routing only; folded/four-pass variants are not yet ported.']
    return dict(distributed_cells=d.rows*d.periods_per_row*nc+(d.rows-1)*d.turn_periods*nc,
                expected_metal_islands=3, turn_max_width_um=float(wm))


def digest(data):
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def generate(design, process, output):
    """Write a GDS plus a self-contained JSON manifest. Does not imply DRC pass.

    The separate validator reads the exported GDS, not these in-memory objects.
    """
    if not isinstance(process,Process) or not isinstance(design,(KITLayout,JosephsonLayout)):
        raise TypeError('Expected a typed layout design and Process')
    design.validate(process)
    c=Canvas(process)
    extra=_kit(c,design) if isinstance(design,KITLayout) else (_josephson(c,design) or {})
    out=Path(output)
    if out.suffix.lower()!='.gds':
        raise ValueError('Output must be a .gds path')
    # Merge only within a layer; preserve layer semantics, including annotations.
    for li in c.ly.layer_indices():
        reg=c.db.Region(c.top.begin_shapes_rec(li)).merged()
        c.top.shapes(li).clear(); c.top.shapes(li).insert(reg)
    die=c.db.Region(c.top.begin_shapes_rec(c.ly.layer(*process.layer('DIE'))))
    for li in c.ly.layer_indices():
        if not (c.db.Region(c.top.begin_shapes_rec(li))-die).is_empty():
            raise ValueError('Generated geometry does not fit inside the die')
    out.parent.mkdir(parents=True,exist_ok=True)
    options=c.db.SaveLayoutOptions()
    options.gds2_max_vertex_count=4000
    c.ly.write(str(out),options)
    pd=process.to_dict(); dd=asdict(design)
    manifest=dict(schema_version=1, generator='twpa-sim.layout.v1',
                  family='kit' if isinstance(design,KITLayout) else design.family,
                  process=pd, process_sha256=digest(pd), design=dd, design_sha256=digest(dd),
                  gds_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
                  components=c.components, pins=c.pins, notes=c.notes, **extra)
    out.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return manifest


def load_design(path, process=None):
    data=_json(Path(path).read_text(encoding='utf-8-sig'))
    _keys(data,('family','geometry'))
    cls=KITLayout if data['family']=='kit' else JosephsonLayout
    if not isinstance(data['geometry'],dict):
        raise ValueError('geometry must be a mapping')
    try:
        if cls is KITLayout:
            return cls.from_process(process,**data['geometry']) if process else cls(**data['geometry'])
        return cls(family=data['family'],**data['geometry'])
    except TypeError as e:
        raise ValueError(f'Unknown/duplicate geometry field: {e}') from e
