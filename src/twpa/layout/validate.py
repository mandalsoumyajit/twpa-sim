"""Independent checks on exported GDS, using KLayout's geometry/DRC engine.

This is research DRC plus an annotated connectivity check, not foundry LVS.
No generator code is imported or replayed to determine polygon connectivity.
"""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from .generate import database, digest
from .process import Process, _json


def validate(gds, manifest=None, report=None, *, drc=True):
    """Check a GDS against its process/design manifest; return a JSON-able report.

    ``report`` is a path stem for .json and .lyrdb outputs. A GDS hash mismatch
    is recorded, permitting inspection of edited layouts against original intent.
    The process/design hashes, unlike the GDS hash, must match the manifest.
    """
    db=database()
    import klayout.rdb as rdb
    gds=Path(gds)
    manifest=Path(manifest) if manifest else gds.with_suffix('.manifest.json')
    m=_json(manifest.read_text(encoding='utf-8-sig'))
    if m['schema_version']!=1 or m['generator']!='twpa-sim.layout.v1':
        raise ValueError('Unsupported layout manifest')
    p=Process(m['process']); pd=p.to_dict()
    if digest(pd)!=m['process_sha256'] or digest(m['design'])!=m['design_sha256']:
        raise ValueError('Manifest process/design hash mismatch')
    ly=db.Layout(); ly.read(str(gds))
    result=dict(schema_version=1, gds_sha256=hashlib.sha256(gds.read_bytes()).hexdigest(),
                process_sha256=m['process_sha256'], design_sha256=m['design_sha256'],
                process_id=p.id, qualification=pd['status'], klayout_version=db.__version__,
                drc_executed=drc, violations=[], rules={}, layers={}, connectivity={})
    result['gds_matches_manifest']=result['gds_sha256']==m['gds_sha256']
    markers=rdb.ReportDatabase('twpa-sim research geometry and DRC')
    rc=markers.create_cell('TWPA')

    def fail(id, message, shapes=None):
        n=shapes.count() if shapes is not None else 1
        if n==0:
            return
        result['violations'].append(dict(id=id,count=n,message=message))
        cat=markers.create_category(id); cat.description=message
        if shapes is not None:
            markers.create_items(rc.rdb_id(),cat.rdb_id(),db.CplxTrans(ly.dbu),shapes)
        else:
            item=markers.create_item(rc.rdb_id(),cat.rdb_id()); item.add_value(message)

    if abs(ly.dbu-p.dbu)>1e-12:
        fail('GEO.DBU','Database unit differs from process')
    if len(ly.top_cells())!=1:
        fail('GEO.TOP','Exactly one top cell is required')
    if not ly.top_cells():
        raise ValueError('GDS has no top cell')
    top=ly.top_cells()[0]
    allowed={tuple(v['gds']) for v in pd['layers'].values()}
    for li in ly.layer_indices():
        info=ly.get_info(li)
        if (info.layer,info.datatype) not in allowed and not db.Region(top.begin_shapes_rec(li)).is_empty():
            fail('GEO.UNKNOWN_LAYER',f'Unexpected layer {info.layer}/{info.datatype}')
        iterator=top.begin_shapes_rec(li)
        while not iterator.at_end():
            shape=iterator.shape()
            if not (shape.is_box() or shape.is_path() or shape.is_polygon() or shape.is_simple_polygon()):
                fail('GEO.UNSUPPORTED_SHAPE','Only polygon/box/path geometry is supported')
                break
            iterator.next()
    raw={name:db.Region(top.begin_shapes_rec(ly.layer(*v['gds']))) for name,v in pd['layers'].items()}
    regs={name:r.merged() for name,r in raw.items()}
    for name,reg in regs.items():
        result['layers'][name]=dict(islands=reg.count(),area_um2=reg.area()*ly.dbu**2)
    if regs['DIE'].count()!=1 or not next(regs['DIE'].each(),db.Polygon()).is_box():
        fail('GEO.DIE','Exactly one rectangular die boundary is required')
    for role,reg in regs.items():
        if role!='DIE':
            fail('GEO.OUTSIDE_DIE',f'{role} outside die',reg-regs['DIE'])
    for role in ('M1',):
        if regs[role].is_empty():
            fail('GEO.MISSING_METAL',f'{role} is empty')
    if drc:
        for rule in pd['rules']:
            name=rule['id']; kind=rule['kind']; ls=rule['layers']
            a=regs[ls[0]]; ticks=round(rule['value_um']/ly.dbu)
            missing_count=0
            if kind=='width':
                bad=a.width_check(ticks,False,db.Region.Euclidian,rule.get('ignore_angle_deg',90))
            elif kind=='space':
                bad=a.space_check(ticks,False,db.Region.Euclidian,rule.get('ignore_angle_deg',90))
            elif kind=='grid':
                r=raw[ls[0]].dup(); r.merged_semantics=False
                bad=r.grid_check(ticks,ticks)
            elif kind=='inside':
                bad=a-regs[ls[1]]
            elif kind=='enclosure':
                # Enclosure distance alone does not detect completely missing outer metal.
                b=regs[ls[1]]
                missing=a-b
                missing_count=missing.count()
                fail(name+'.MISSING','Inner geometry is not contained in outer layer',missing)
                bad=a.enclosed_check(b,ticks)
            else:
                raise ValueError(f'Unimplemented rule {kind}')
            result['rules'][name]=dict(kind=kind,status=rule['status'],violations=bad.count()+missing_count)
            if 'ignore_angle_deg' in rule:
                result['rules'][name]['ignore_angle_deg']=rule['ignore_angle_deg']
            fail(name,f'{kind}: {ls}, limit {rule["value_um"]} um',bad)
    if m['family']=='kit':
        # Signal + two ground combs, with all three net assignments verified below.
        if regs['M1'].count()!=3:
            fail('GEO.KIT_ISLANDS',f'Expected 3 conductor islands, found {regs["M1"].count()}')
        if regs['SIG_PAD'].count()!=2:
            fail('GEO.KIT_PORTS','Expected two signal landing pads')
        for role in ('SIG_PAD','GND_PAD'):
            fail('GEO.PAD_METAL',f'{role} outside NbN',regs[role]-regs['M1'])
        if not (regs['SIG_PAD'] & regs['GND_PAD']).is_empty():
            fail('GEO.PAD_SHORT','Ground/signal pad masks overlap')
    else:
        expected=Counter()
        d=m['design']; n=d['cells']; family=m['family']
        expected['JJ']=n*(1+(d['large_junctions'] if family=='snail3wm' else 1 if family=='squid_jpa' else 0))
        expected['CAP']=1 if family in ('jpa','squid_jpa') else n+1
        expected['IND']=n if family=='rfsquid3wm' else 0
        for kind,count in expected.items():
            if regs[kind].count()!=count or sum(c['kind']==kind for c in m['components'])!=count:
                fail('GEO.COMPONENT_COUNT',f'Expected {count} {kind} components')
        # Every annotated component must match an actual polygon, not only a count.
        for kind in ('JJ','CAP','IND'):
            intended=db.Region()
            for c in m['components']:
                if c['kind']==kind:
                    intended.insert(db.Box(*(round(v/ly.dbu) for v in c['box_um'])))
            fail('GEO.COMPONENT_SHAPE',f'{kind} differs from circuit intent',regs[kind]^intended)
        fail('GEO.BARRIER_CONTACT','Galvanic contacts overlap tunnel/capacitor windows',
             regs['VIA'] & (regs['JJ']+regs['CAP']))

    # Cut annotated inductors only in the extraction view, never in exported metal.
    extraction={'M1':regs['M1']-regs.get('IND',db.Region())}
    if 'M2' in regs:
        extraction['M2']=regs['M2']
    pieces={}; parent={}
    for role,reg in extraction.items():
        for i,poly in enumerate(reg.merged().each()):
            key=(role,i); pieces[key]=db.Region(poly); parent[key]=key

    def root(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]]; x=parent[x]
        return x

    def hits(region,role):
        return [key for key,poly in pieces.items() if key[0]==role and not (poly & region).is_empty()]

    if 'VIA' in regs:
        for via in regs['VIA'].each():
            vr=db.Region(via)
            low=hits(vr,'M1'); high=hits(vr,'M2')
            if len(low)!=1 or len(high)!=1:
                fail('GEO.VIA','Contact must connect exactly one M1 and one M2 island')
            else:
                parent[root(high[0])]=root(low[0])
    assignments=[]
    for pin in m['pins']:
        x,y=(round(v/ly.dbu) for v in pin['at_um'])
        found=hits(db.Region(db.Box(x-1,y-1,x+1,y+1)),pin['layer'])
        if len(found)!=1:
            fail('GEO.PIN',f'Pin {pin["name"]} is absent or ambiguous')
        else:
            assignments.append((pin['net'],root(found[0])))
    for c in m['components']:
        box=db.Box(*(round(v/ly.dbu) for v in c['box_um'])); cr=db.Region(box)
        if c['kind']=='IND':
            found=hits(cr.sized(1),'M1')
            found.sort(key=lambda k:pieces[k].bbox().center().x)
            segment=(cr & regs['M1']).merged()
            if (len(found)!=2 or segment.count()!=1 or
                    segment.bbox().left!=box.left or segment.bbox().right!=box.right):
                fail('GEO.INDUCTOR','Inductor annotation must bridge two cut M1 islands')
                continue
            terminals=found
        else:
            low=hits(cr,'M1'); high=hits(cr,'M2')
            if len(low)!=1 or len(high)!=1:
                fail('GEO.COMPONENT_TERMINALS',f'{c["id"]} must join one lower and one upper electrode')
                continue
            terminals=low+high if c['kind']=='JJ' else high+low
        for net,key in zip(c['nets'],terminals):
            assignments.append((net,root(key)))
    nets=defaultdict(set); physical=defaultdict(set)
    for net,key in assignments:
        nets[net].add(key); physical[key].add(net)
    for net,keys in nets.items():
        if len(keys)!=1:
            fail('GEO.OPEN',f'Intended net {net} spans {len(keys)} disconnected metal nets')
    for key,names in physical.items():
        if len(names)!=1:
            fail('GEO.SHORT',f'Physical metal shorts intended nets {sorted(names)}')
    floating={root(key) for key in pieces}-set(physical)
    if floating:
        fail('GEO.FLOATING',f'{len(floating)} unassigned conductor islands')
    result['connectivity']=dict(physical_nets=len({root(key) for key in pieces}),
                                intended_nets=len(nets),unassigned_islands=len(floating))
    result['passed']=not result['violations']
    result['scope']='Research geometry and explicit rule subset; no foundry signoff, EM extraction, or full LVS.'
    if report:
        stem=Path(report); stem.parent.mkdir(parents=True,exist_ok=True)
        stem.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        markers.save(str(stem.with_suffix('.lyrdb')))
    return result
