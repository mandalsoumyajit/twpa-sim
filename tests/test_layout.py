"""Export/reimport, independent geometry, real DRC and deliberate defects."""
from copy import deepcopy
from pathlib import Path
import os
import tempfile
import unittest

from twpa.layout import Process, load_process

try:
    import klayout.db as db
except ImportError:
    db=None


class ProcessTests(unittest.TestCase):
    def test_profiles_and_provenance(self):
        kit=load_process(); jj=load_process('jj_overlap_illustrative')
        self.assertEqual(kit.layer('M1'),(2,0))
        self.assertEqual(jj.to_dict()['status'],'illustrative')
        self.assertEqual(kit.parameter('sheet_inductance','pH/sq'),20)
        self.assertTrue(all(r['status']=='engineering_assumption' for r in kit.to_dict()['rules']))
        self.assertEqual(kit.coordinate(-.58),-580)
        with self.assertRaises(ValueError): kit.coordinate(.5801)
        with self.assertRaises(ValueError): kit.parameter('sheet_inductance','H')
        d=kit.to_dict(); d['layers']['M1']['gds'][0]=9
        self.assertEqual(kit.layer('M1'),(2,0))

    def test_invalid_contracts(self):
        original=load_process().to_dict()
        def damage(d,case):
            if case==0: d['unexpected']=1
            if case==1: d['layers']['M1']['gds'][1]=None
            if case==2: d['layers']['M1']['gds']=d['layers']['DIE']['gds']
            if case==3: d['rules'][0]['layers']=['unknown']
            if case==4: d['rules'][0]['kind']='not_implemented'
            if case==5: d['rules'][0]['value_um']=float('nan')
            if case==6: d['grid_um']=.0015
            if case==7: d['rules'][0]['source']='no_source'
            if case==8: d['status']='signoff'
            if case==9: d['rules'][0]['value_um']=.2001
            if case==10: d['layers']['M1']['gds'][0]=True
            if case==11: d['rules'][0]['ignore_angle_deg']=float('nan')
            if case==12: d['parameters']['pitch']['unit']=''
            if case==13: d['rules'].append(deepcopy(d['rules'][0]))
        for case in range(14):
            with self.subTest(case=case):
                d=deepcopy(original); damage(d,case)
                with self.assertRaises(ValueError): Process(d)


@unittest.skipIf(db is None,'Install twpa-sim[layout]')
class LayoutTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.jj=load_process('jj_overlap_illustrative')

    def tearDown(self):
        self.temp.cleanup()

    def build(self,family='jpa',**kwargs):
        from twpa.layout.generate import generate, JosephsonLayout, KITLayout
        out=self.root/(family+'.gds')
        design=KITLayout.from_process(load_process(),rows=2,periods_per_row=2,**kwargs) if family=='kit' else JosephsonLayout(family=family,**kwargs)
        m=generate(design,load_process() if family=='kit' else self.jj,out)
        return out,m

    def check(self,path):
        from twpa.layout.validate import validate
        return validate(path,report=self.root/(path.stem+'-check'))

    def mutate(self,path,action):
        ly=db.Layout(); ly.read(str(path)); top=ly.top_cell()
        action(ly,top)
        ly.write(str(path))

    def test_all_families_geometry_and_drc(self):
        for family in ('kit','jpa','squid_jpa','jtwpa4wm','rfsquid3wm','snail3wm'):
            with self.subTest(family=family):
                path,m=self.build(family,**({} if family=='kit' or 'jpa' in family else {'cells':3}))
                r=self.check(path)
                self.assertTrue(r['passed'],r['violations'])
                self.assertTrue(r['gds_matches_manifest'])
                self.assertEqual(r['connectivity']['unassigned_islands'],0)
                if family=='kit':
                    self.assertEqual(m['distributed_cells'],378)
                else:
                    jj=next(c for c in m['components'] if c['kind']=='JJ')
                    area=jj['estimates']['drawn_area_um2']
                    self.assertAlmostEqual(jj['estimates']['Ic_A'],area*5e-6)
                    self.assertAlmostEqual(jj['estimates']['Cj_F'],area*60e-15)
                    if 'jpa' not in family:
                        self.assertAlmostEqual(sum(c['estimates']['C_F'] for c in m['components'] if c['kind']=='CAP'),
                                               m['design']['cells']*144e-15,delta=1e-26)

    def test_changed_geometry(self):
        from twpa.layout.generate import KITLayout,generate
        path=self.root/'changed.gds'
        generate(KITLayout(rows=3,periods_per_row=3,spine_width=.64,loaded_cells=8,
                           finger_length=65,turn_radius=160),load_process(),path)
        self.assertTrue(self.check(path)['passed'])
        path,_=self.build('snail3wm',cells=2,large_junctions=4,junction_width=2.5,
                          small_junction_width=.8,capacitor_width=10,wire_width=1.5)
        r=self.check(path)
        self.assertTrue(r['passed'],r['violations'])

    def test_original_kit_primitive_fixture(self):
        import json
        from twpa.layout.generate import KITLayout
        from twpa.layout.kit import primitives
        ref=json.loads((Path(__file__).parents[1]/'validation/layout-reference-primitives.json').read_text())
        ly=db.Layout(); ly.dbu=.001
        built=primitives(ly,KITLayout())
        for name,polys in zip(('straight','turn_right','turn_left','launch_left','launch_right'),built[:5]):
            self.assertAlmostEqual(db.Region(polys).area()*1e-6,ref['primitives'][name]['area_um2'],delta=1e-6)
            self.assertEqual(db.Region(polys).merged().count(),ref['primitives'][name]['islands'])
        self.assertAlmostEqual(built[-1],ref['turn_max_width_um'])

    def test_electrical_sizing_and_incompatible_capacitance(self):
        from twpa.layout.bindings import from_circuit
        from twpa.jpa_circuit import ParallelJPA
        from twpa.jtwpa import Chain,ResonantLoad
        p=self.jj
        jpa=ParallelJPA(20e-6,144e-15,240e-15)
        d=from_circuit(jpa,p)
        self.assertAlmostEqual(d.junction_width,2)
        self.assertAlmostEqual(d.capacitor_height,12)
        # Larger resonant shunts are supported without changing the junction area.
        large=from_circuit(ParallelJPA(20e-6,12e-12,240e-15),p)
        self.assertAlmostEqual(large.capacitor_height,1000)
        chain=Chain(20e-6,144e-15,240e-15,8,cell_length=60e-6)
        d=from_circuit(chain,p,cell_pitch=60)
        self.assertEqual(d.cells,8)
        with self.assertRaisesRegex(ValueError,'incompatible'):
            from_circuit(ParallelJPA(20e-6,144e-15,1e-12),p)
        with self.assertRaisesRegex(ValueError,'additional layout primitives'):
            from_circuit(Chain(20e-6,144e-15,240e-15,8,resonator=ResonantLoad(1e-15,1e-12,1e-9)),p)

    def test_shorts_and_opens(self):
        path,_=self.build()
        def short(ly,top):
            top.shapes(ly.layer(10,0)).insert(db.Box(21000,-1000,35000,1000))
        self.mutate(path,short)
        r=self.check(path)
        self.assertIn('GEO.SHORT',{v['id'] for v in r['violations']})
        self.assertFalse(r['gds_matches_manifest'])
        path,_=self.build()
        def cut(ly,top):
            layer=ly.layer(10,0); region=db.Region(top.begin_shapes_rec(layer))
            region-=db.Region(db.Box(10000,-4000,11000,4000))
            top.shapes(layer).clear(); top.shapes(layer).insert(region)
        self.mutate(path,cut)
        self.assertIn('GEO.OPEN',{v['id'] for v in self.check(path)['violations']})

    def test_component_contact_and_inductor_defects(self):
        for role,code in ((40,'GEO.COMPONENT_COUNT'),(30,'GEO.OPEN')):
            path,_=self.build()
            self.mutate(path,lambda ly,top: top.shapes(ly.layer(role,0)).clear())
            r=self.check(path)
            self.assertIn(code,{v['id'] for v in r['violations']})
        path,_=self.build('rfsquid3wm',cells=2)
        self.mutate(path,lambda ly,top: top.shapes(ly.layer(60,0)).clear())
        self.assertFalse(self.check(path)['passed'])
        path,_=self.build('rfsquid3wm',cells=2)
        def cut_inside_inductor(ly,top):
            layer=ly.layer(10,0)
            r=db.Region(top.begin_shapes_rec(layer))-db.Region(db.Box(29500,17000,30500,23000))
            top.shapes(layer).clear(); top.shapes(layer).insert(r)
        self.mutate(path,cut_inside_inductor)
        self.assertIn('GEO.INDUCTOR',{v['id'] for v in self.check(path)['violations']})

    def test_drc_defects(self):
        cases=[('width',db.Box(10000,6000,10100,9000),'M1.W',10),
               ('space',db.Box(10000,1100,12000,3100),'M1.S',10),
               ('enclosure',db.Box(15000,8000,16000,9000),'VIA.M1.MISSING',30),
               ('die',db.Box(-100000,-100000,-90000,-90000),'GEO.OUTSIDE_DIE',10),
               ('layer',db.Box(0,0,1000,1000),'GEO.UNKNOWN_LAYER',99)]
        for name,box,code,layer in cases:
            with self.subTest(defect=name):
                path,_=self.build()
                self.mutate(path,lambda ly,top:top.shapes(ly.layer(layer,0)).insert(box))
                self.assertIn(code,{v['id'] for v in self.check(path)['violations']})

    def test_offgrid_and_malformed_manifest(self):
        from twpa.layout.generate import generate,JosephsonLayout
        from twpa.layout.validate import validate
        import json
        d=self.jj.to_dict(); d['grid_um']=.005
        for r in d['rules']:
            if r['kind']=='grid': r['value_um']=.005
        path=self.root/'grid.gds'
        generate(JosephsonLayout(),Process(d),path)
        self.mutate(path,lambda ly,top:top.shapes(ly.layer(10,0)).insert(db.Box(10001,6000,12001,8000)))
        self.assertIn('M1.GRID',{v['id'] for v in self.check(path)['violations']})
        data=json.loads(path.with_suffix('.manifest.json').read_text())
        data['process']['rules'][0]['value_um']=.010
        path.with_suffix('.manifest.json').write_text(json.dumps(data))
        with self.assertRaises(ValueError): validate(path)

    @unittest.skipUnless(os.environ.get('TWPA_KLAYOUT_EXECUTABLE'),'Optional native Ruby DRC engine')
    def test_native_ruby_drc(self):
        from twpa.layout.native_drc import run_klayout_drc
        executable=os.environ['TWPA_KLAYOUT_EXECUTABLE']
        for family in ('kit','jpa','squid_jpa','jtwpa4wm','rfsquid3wm','snail3wm'):
            path,_=self.build(family)
            r=run_klayout_drc(path,load_process() if family=='kit' else self.jj,
                              self.root/(family+'.lyrdb'),executable)
            self.assertTrue(r['passed'],r)
        path,_=self.build()
        defects=[(10,db.Box(10000,6000,10100,9000)),
                 (10,db.Box(10000,1100,12000,3100)),
                 (30,db.Box(15000,8000,16000,9000))]
        for i,(layer,box) in enumerate(defects):
            path,_=self.build()
            self.mutate(path,lambda ly,top:top.shapes(ly.layer(layer,0)).insert(box))
            self.assertFalse(run_klayout_drc(path,self.jj,self.root/f'bad{i}.lyrdb',executable)['passed'])
        from twpa.layout.generate import generate,JosephsonLayout
        data=self.jj.to_dict(); data['grid_um']=.005
        for r in data['rules']:
            if r['kind']=='grid': r['value_um']=.005
        process=Process(data)
        generate(JosephsonLayout(),process,path)
        self.mutate(path,lambda ly,top:top.shapes(ly.layer(10,0)).insert(db.Box(10001,6000,12001,8000)))
        self.assertFalse(run_klayout_drc(path,process,self.root/'bad-grid.lyrdb',executable)['passed'])


if __name__=='__main__':
    unittest.main()
