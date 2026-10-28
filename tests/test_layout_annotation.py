"""Finite-network physics, geometry binding, and calibrated KIT reduction."""
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
import numpy as np
from twpa.em import LayoutNetwork,load_circuit,kit_calibration,back_annotate_kit
from twpa.jpa_circuit import PHIBAR,Tone,Transient
from twpa.models import Bias


def rf_network():
    # A physical loop: one linear arm and one RCSJ arm between x and ground.
    return dict(family='rf-test',nodes=['x','g'],reference='g',
        metal_branches=[dict(positive='x',negative='g')],inductance_H=[[80e-12]],
        junctions=[dict(id='J',positive='x',negative='g',critical_current=2e-6,capacitance=60e-15,resistance=1e6)],
        capacitors=[dict(id='C',positive='x',negative='g',capacitance=150e-15,resistance=1e6)],
        ports=[dict(positive='x',negative='g',impedance=50.)])


class TerminalNetworkPhysics(unittest.TestCase):
    def test_total_l_certificate_keeps_fixed_gates(self):
        from twpa.em import certify_magnetic
        from twpa.layout.generate import digest
        def sample(**settings):
            s=dict(pitch_um=2.,local_pitch_um=.125,nh=2,patch_side_um=1.,lambda_scale=1.,frequencies_hz=[1e6]);s.update(settings)
            d=dict(schema_version=1,kind='layout_london_network',source={'fixture':'SYNTHETIC'},manifest_sha256='test',
                stack={},patches=[],reduction=[[1.,0.],[0.,1.]],settings=s,
                ports=[dict(name='a',positive_patch='A',negative_patch='R'),dict(name='b',positive_patch='B',negative_patch='R')],
                run=dict(simulator=dict(binary_sha256='SYNTHETIC')),
                samples=[dict(frequency_Hz=s['frequencies_hz'][0],L_H=[[100e-12,10e-12],[10e-12,100e-12]])])
            d['artifact_sha256']=digest(d);return d
        p=sample();c=dict(bulk=sample(pitch_um=3.),thickness=sample(nh=1),contact=sample(local_pitch_um=.25),frequency=sample(frequencies_hz=[1e11]))
        accepted=certify_magnetic(p,c)
        self.assertTrue(accepted['certification']['passed']);self.assertNotIn('kinetic_energy_H',accepted)
        bad=deepcopy(c);d=bad['contact'];d['samples'][0]['L_H'][0][1]=d['samples'][0]['L_H'][1][0]=10.6e-12
        d.pop('artifact_sha256');d['artifact_sha256']=digest(d)
        rejected=certify_magnetic(p,bad)
        self.assertFalse(rejected['certification']['passed'])
        self.assertEqual(rejected['certification']['absolute_coupling_gate'],.005)

    def test_rf_screening_and_flux_orientation(self):
        c=LayoutNetwork(rf_network());bias=.12;theta=c.equilibrium(bias)[0]
        # Closing-JJ orientation defines the signed loop flux.
        beta=c.L[0,0]*c.ic[0]/PHIBAR
        self.assertAlmostEqual(theta+beta*np.sin(theta),2*np.pi*bias,places=11)
        np.testing.assert_allclose(np.hstack((c.BL,c.BJ))@c.cycles.T,0,atol=1e-14)
        np.testing.assert_allclose(c.CL@c.flux_vector(bias),-2*np.pi*PHIBAR*bias)

    def test_zero_pump_matches_independent_linear_impedance(self):
        c=LayoutNetwork(rf_network());f=5e9;w=2*np.pi*f
        Y=1/c.R[0]+2e-6+1j/(w*c.L[0,0])+1j*c.ic[0]/(w*PHIBAR)-1j*w*210e-15
        expected=2/(50*Y)-1
        self.assertAlmostEqual(c.passive([f])[0,0],expected,places=12)
        p=c.periodic(6e9,0.);r=c.response(p,[f],sidebands=2)
        self.assertAlmostEqual(r['waves'][0,2,0],expected,places=11)
        self.assertLess(np.max(abs(r['waves'][0,[0,1,3,4],:])),1e-12)
        self.assertLessEqual(abs(expected)**2,1.)

    def test_signed_orientation_invariance_and_input_ownership(self):
        d=rf_network();c=LayoutNetwork(d);d['junctions'][0]['critical_current']=1.
        self.assertEqual(c.ic[0],2e-6)
        d=rf_network();d['metal_branches'][0]=dict(positive='g',negative='x')
        b=LayoutNetwork(d)
        np.testing.assert_allclose(c.passive([4e9,7e9]),b.passive([4e9,7e9]),rtol=1e-13)

    def test_deck_keeps_junction_capacitance_once_and_signed_mutuals(self):
        d=rf_network();d['nodes'].append('y');d['metal_branches'].append(dict(positive='y',negative='g'))
        d['inductance_H']=[[80e-12,-4e-12],[-4e-12,40e-12]]
        d['capacitors'].append(dict(id='C2',positive='y',negative='g',capacitance=100e-15,resistance=1e6))
        c=LayoutNetwork(d);deck=c.netlist([Tone(5e9,1e-12)],Transient(1e-14,1e-9,0.,1e-12),flux_bias=.05)
        self.assertIn('KM1_0 LM1 LM0 -',deck);self.assertEqual(deck.count('.model JJ'),1)
        self.assertEqual(deck.count('CC'),2);self.assertIn('IFLUX0',deck)
        basis=c.flux_vector([1.])[:,None]/1e-3
        coil=max(1e-9,2*np.linalg.eigvalsh(basis.T@c.invL@basis)[-1])
        self.assertGreater(np.linalg.eigvalsh(np.block([[c.L,basis],[basis.T,np.array([[coil]])]]))[0],0)

    def test_reject_invalid_or_unresolved_physics(self):
        for key,value in [('inductance_H',[[-1.]]),('reference','absent')]:
            d=rf_network();d[key]=value
            with self.assertRaises(ValueError):LayoutNetwork(d)
        c=LayoutNetwork(rf_network())
        with self.assertRaises(ValueError):c.flux_vector([0.,0.])
        with self.assertRaises(ValueError):c.response(c.periodic(5e9,0.),[5e9])
        with self.assertRaises(ValueError):c.response(c.periodic(5e9,0.),[4e9],quadrature=1)


try:
    import klayout.db
    HAVE_LAYOUT=True
except ImportError:HAVE_LAYOUT=False


@unittest.skipUnless(HAVE_LAYOUT,'Install the layout extra')
class LayoutAnnotationGeometry(unittest.TestCase):
    def test_complete_annotation_binding_and_reuse_contracts(self):
        # Deliberately synthetic solver values exercise contracts only. They
        # are never saved as native validation evidence.
        import json
        from twpa.layout import load_process
        from twpa.layout.generate import generate,JosephsonLayout,digest
        from twpa.em import layout_magnetic_job,capacitor_job,load_stack,back_annotate_layout,reuse_capacitor,read_extraction
        def seal(d):
            d=deepcopy(d);d.pop('artifact_sha256',None);d['artifact_sha256']=digest(d);return d
        with tempfile.TemporaryDirectory() as tmp:
            for family,cells,loops in [('squid_jpa',1,1),('jtwpa4wm',2,0),('rfsquid3wm',2,2),('snail3wm',1,1)]:
                path=Path(tmp)/(family+'.gds');m=generate(JosephsonLayout(family=family,cells=cells),load_process('jj_overlap_illustrative'),path)
                job=layout_magnetic_job(path,pitch_um=4,local_pitch_um=.5,nh=1,max_filaments=40000,frequencies_hz=(1e6,))
                field={k:deepcopy(v) for k,v in job.items() if k not in ('deck','segments','coordinates')}
                field.update(run=dict(job_sha256=digest(job)),inductance_H=(np.eye(len(job['ports']))*1e-11).tolist(),certification=dict(passed=True))
                field=seal(field);electric={}
                for c in m['components']:
                    if c['kind']!='CAP':continue
                    e=capacitor_job(path,load_stack(),component_id=c['id'])
                    cap=(c['box_um'][3]-c['box_um'][1])/12*100e-15
                    e.update(simulator=dict(name='openEMS',version='SYNTHETIC'),fit_accepted=True,
                             capacitance_F=cap,capacitance_by_frequency_F=[cap]*6,convergence=dict(passed=True))
                    electric[c['id']]=seal(e)
                output=Path(tmp)/(family+'.json')
                circuit,record=back_annotate_layout(path,electric,field,output=output)
                self.assertEqual(len(circuit.cycles),loops)
                self.assertEqual(len(circuit.ic),sum(c['kind']=='JJ' for c in m['components']))
                self.assertEqual(len(circuit.data['capacitors']),len(electric))
                np.testing.assert_allclose([c['capacitance'] for c in circuit.data['capacitors']],
                                           [d['capacitance_F'] for d in electric.values()])
                np.testing.assert_allclose(load_circuit(output).L,circuit.L)
                with self.assertRaises(ValueError):back_annotate_layout(path,{},field)
                failed=deepcopy(field);failed['certification']['passed']=False
                with self.assertRaises(ValueError):back_annotate_layout(path,electric,seal(failed))
                names=list(electric)
                reused=reuse_capacitor(path,names[0],electric[names[0]])
                self.assertIn('reused_native_extraction',read_extraction(reused))
                reused['capacitance_F']*=2
                with self.assertRaises(ValueError):read_extraction(seal(reused))
                if len(names)>1:
                    wrong={k:electric[names[0]] for k in names}
                    with self.assertRaises(ValueError):back_annotate_layout(path,wrong,field)
                stale=deepcopy(field);stale['source']['gds_sha256']='stale'
                with self.assertRaises(ValueError):back_annotate_layout(path,electric,seal(stale))

    def test_all_josephson_terminals_and_distinct_end_capacitors(self):
        from twpa.layout import load_process
        from twpa.layout.generate import generate,JosephsonLayout
        from twpa.em import layout_magnetic_job,capacitor_job,load_stack
        with tempfile.TemporaryDirectory() as tmp:
            for family in ('squid_jpa','jtwpa4wm','rfsquid3wm','snail3wm'):
                path=Path(tmp)/(family+'.gds')
                generate(JosephsonLayout(family=family),load_process('jj_overlap_illustrative'),path)
                job=layout_magnetic_job(path,pitch_um=4,local_pitch_um=.5,nh=1,max_filaments=40000,frequencies_hz=(1e6,))
                self.assertTrue(job['audit']['physical_ports_connected'])
                self.assertEqual(np.linalg.matrix_rank(job['reduction']),len(job['ports']))
                names={p['name'] for p in job['patches']}
                self.assertTrue({'INPUT','GROUND0'}<=names)
                import json
                m=json.loads(path.with_suffix('.manifest.json').read_text())
                for comp in m['components']:
                    if comp['kind'] in ('JJ','CAP'):
                        self.assertTrue({comp['id']+'_M1',comp['id']+'_M2'}<=names)
                caps=[c for c in m['components'] if c['kind']=='CAP']
                if len(caps)>1:
                    with self.assertRaises(ValueError):capacitor_job(path,load_stack())
                    for cap in caps:
                        e=capacitor_job(path,load_stack(),component_id=cap['id'])
                        self.assertEqual(e['source']['component_id'],cap['id'])

    def test_kit_reference_reduction_and_stale_calibration(self):
        from twpa.layout import load_process
        from twpa.layout.generate import generate,KITLayout
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'kit.gds'
            generate(KITLayout(rows=2,periods_per_row=1),load_process(),path)
            cal=kit_calibration(path,capacitance_per_finger_F=(11.16e-15,4.19e-15),
                stub_inductance_H=(20e-12*(75.5/.5)/3,20e-12*(27.5/.5)/3),temperature_K=.02,
                source='Original research design coupon calibration; example only',turn_reduction='sheet_squares_matched_periodic')
            model,record=back_annotate_kit(path,cal,Bias(),output=Path(tmp)/'circuit.json')
            self.assertEqual(model.design.repeats,5);self.assertEqual(model.design.total_cells,270)
            self.assertAlmostEqual(model.design.sections[0].cell.inductance,20e-12*(1.5/.58)*1.04)
            self.assertLess(record['checks']['turn_sheet_relative_error'],1e-8)
            restored=load_circuit(Path(tmp)/'circuit.json')
            np.testing.assert_allclose(restored.linear(np.linspace(0,4e9,81)).s21,model.linear(np.linspace(0,4e9,81)).s21)
            passive=model.linear(np.linspace(0,4e9,81))
            np.testing.assert_allclose(abs(passive.s11)**2+abs(passive.s21)**2,1.,atol=1e-10)
            from twpa.models import TLS
            grid=np.linspace(0,8e9,801)
            analytic=model.small_signal([2e9],frequency=grid,pump_frequency=6e9)
            coupled=model.coupled([2e9],TLS(),frequency=grid,pump_frequency=6e9)
            np.testing.assert_allclose(coupled.gain,analytic.gain,rtol=1e-8)
            np.testing.assert_allclose(coupled.added_noise,(1-1/coupled.gain)/2,atol=1e-10)
            with self.assertRaises(ValueError):back_annotate_kit(path,cal,Bias(dc_current=.1e-3))
            damaged=deepcopy(cal);damaged['capacitance_per_finger_F'][0]*=2
            with self.assertRaises(ValueError):back_annotate_kit(path,damaged,Bias())


if __name__=='__main__':unittest.main()
