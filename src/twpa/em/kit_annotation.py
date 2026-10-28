"""Checked KIT geometry to the existing calibrated periodic propagation model.

This is a process/calibration reduction, not a full-serpentine field extraction.
Turns use the original generator's biased sheet-squares matching assumption;
capacitance, stub response, participation and material temperature are explicit.
"""
from dataclasses import dataclass,asdict
import hashlib
import json
from pathlib import Path
import numpy as np
from ..models import Cell,Section,Design,Bias,Nonlinearity
from ..layout.generate import digest
from ..layout.process import _json,Process,positive
from ..layout.validate import validate


def _layout(gds):
    p=Path(gds);mp=p.with_suffix('.manifest.json');m=_json(mp.read_text(encoding='utf-8'))
    if m['family']!='kit':raise ValueError('KIT layout required')
    checked=validate(p)
    if not checked['passed'] or not checked['gds_matches_manifest']:raise ValueError('Checked matching KIT GDS required')
    source=dict(gds_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),manifest_sha256=hashlib.sha256(mp.read_bytes()).hexdigest(),
                process_sha256=m['process_sha256'],design_sha256=m['design_sha256'])
    return m,source


def kit_calibration(gds,*,capacitance_per_finger_F,temperature_K,source,
                    stub_inductance_H=(0.,0.),participation=(0.,0.),sheet_inductance_H_per_square=None,
                    turn_reduction):
    """Bind explicit unloaded/loaded calibration pairs to a particular layout.

    The required turn policy declares the approximation; no capacitance is
    inferred by multiplying an old calibration by a changed finger length.
    """
    m,binding=_layout(gds)
    if turn_reduction!='sheet_squares_matched_periodic':raise ValueError('Explicit supported turn-reduction policy required')
    if not isinstance(source,str) or not source.strip():raise ValueError('Calibration source description required')
    positive(temperature_K,'calibration temperature')
    values={}
    for name,items in [('capacitance_per_finger_F',capacitance_per_finger_F),('stub_inductance_H',stub_inductance_H),('participation',participation)]:
        a=np.asarray(items,float)
        if a.shape!=(2,) or np.any(~np.isfinite(a)) or np.any(a<0) or (name=='capacitance_per_finger_F' and np.any(a<=0)):
            raise ValueError('Finite unloaded/loaded calibration pairs required')
        if name=='participation' and np.any(a>1):raise ValueError('Participation must lie in [0,1]')
        values[name]=a.tolist()
    sheet=sheet_inductance_H_per_square
    if sheet is None:sheet=Process(m['process']).parameter('sheet_inductance','pH/sq')*1e-12
    positive(sheet,'zero-bias sheet inductance')
    d=dict(schema_version=1,kind='kit_layout_calibration',status='user_supplied_not_EM_certificate',binding=binding,
           temperature_K=temperature_K,source_description=source,sheet_inductance_H_per_square=sheet,
           turn_reduction=turn_reduction,**values)
    d['artifact_sha256']=digest(d);return d


@dataclass(frozen=True)
class KITPeriodicModel:
    design: Design
    bias: Bias

    def linear(self,frequency):
        from ..linear import linear
        return linear(self.design,frequency,self.bias.impedance)

    def small_signal(self,signal_frequency,**settings):
        from ..amplifier import small_signal
        return small_signal(self.design,self.bias,signal_frequency,**settings)

    def coupled(self,signal_frequency,tls,**settings):
        from ..amplifier import coupled
        return coupled(self.design,self.bias,signal_frequency,tls,**settings)


def back_annotate_kit(gds,calibration,bias,*,output=None):
    m,binding=_layout(gds)
    d=dict(calibration) if isinstance(calibration,dict) else _json(Path(calibration).read_text(encoding='utf-8'))
    h=d.pop('artifact_sha256',None)
    if h!=digest(d) or d.get('kind')!='kit_layout_calibration':raise ValueError('Calibration hash/type mismatch')
    if d['binding']!=binding:raise ValueError('Stale KIT calibration; bind new geometry explicitly')
    if not isinstance(bias,Bias):raise TypeError('Explicit electrical Bias required')
    g=m['design'];ratio=bias.dc_current/bias.scaling_current
    if not np.isclose(ratio,g['bias_ratio'],rtol=0,atol=1e-10):raise ValueError('Electrical DC bias differs from geometry turn-matching bias')
    sheet=d['sheet_inductance_H_per_square'];lc=sheet*g['pitch']/g['spine_width']*(1+ratio**2)
    cells=[Cell(lc,c,ls,p) for c,ls,p in zip(d['capacitance_per_finger_F'],d['stub_inductance_H'],d['participation'])]
    repeats=g['rows']*g['periods_per_row']+(g['rows']-1)*g['turn_periods']
    design=Design((Section(cells[0],g['unloaded_cells']//2),Section(cells[1],g['loaded_cells']),
                   Section(cells[0],g['unloaded_cells']//2)),repeats,name='layout-'+m['process']['id'])
    if design.total_cells!=m['distributed_cells']:raise ValueError('Layout/model distributed cell count differs')
    count=(g['unloaded_cells']+g['loaded_cells'])*g['turn_periods'];theta=np.pi/2-(np.arange(count)+.5)*np.pi/count
    distance=np.minimum(np.pi/2-theta,theta+np.pi/2);u=np.clip(distance/(g['taper_cells']*np.pi/count),0,1)
    width=g['spine_width']+(m['turn_max_width_um']-g['spine_width'])*u*u*(3-2*u)
    turn_l=sheet*(np.pi*g['turn_radius']/count)/width*(1+ratio**2*(g['spine_width']/width)**2)
    error=float(abs(np.sum(turn_l)/(count*lc)-1))
    if error>1e-8:raise ValueError('Turn sheet-squares matching failed')
    model=KITPeriodicModel(design,bias)
    record=dict(schema_version=1,circuit_type='KITPeriodicModel',circuit=dict(design=asdict(design),bias=asdict(bias)),
                source=binding,calibration_sha256=h,calibration=d,
                checks=dict(geometry_passed=True,cell_count=design.total_cells,turn_sheet_relative_error=error),
                route=dict(straight_cells=g['rows']*g['periods_per_row']*(g['unloaded_cells']+g['loaded_cells']),
                           turn_cells=(g['rows']-1)*count,turn_inductances_at_bias_H=turn_l.tolist()),
                scope='Process/calibration-backed periodic KIT reduction, not a full-layout EM extraction.',
                limitations=['Straight-cell capacitance/stub/participation are explicitly supplied calibration values; turns use their declared periodic reuse approximation.',
                    'Biased sheet-squares matching preserves total turn inductance but not detailed bend dispersion or local nonlinear-current variation.',
                    'Launches, separate ground-comb package connections, radiation and nonlocal mutual inductance remain external to this reduction.'])
    record['artifact_sha256']=digest(record)
    if output:
        dest=Path(output);dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return model,record


def load_kit_data(data):
    d=data['design'];b=dict(data['bias'])
    if b.get('nonlinearity') is not None:b['nonlinearity']=Nonlinearity(**b['nonlinearity'])
    def sections(values):return tuple(Section(Cell(**s['cell']),s['count']) for s in values)
    return KITPeriodicModel(Design(sections(d['sections']),d['repeats'],sections(d['tail']),d['name']),Bias(**b))
