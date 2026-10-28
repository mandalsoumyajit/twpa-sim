"""Independent native-field and JoSIM fixtures for layout-derived networks."""
import hashlib
import json
from pathlib import Path
import unittest
import numpy as np
from twpa.em import LayoutNetwork
from twpa.layout.generate import digest


class LayoutNetworkNativeReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((Path(__file__).resolve().parents[1]/'validation/layout-network-reference.json').read_text(encoding='utf-8'))

    def test_complete_native_evidence_and_signed_matrices(self):
        data=dict(self.data);seal=data.pop('artifact_sha256');self.assertEqual(digest(data),seal)
        self.assertEqual({c['family'] for c in data['cases'].values()},{'squid_jpa','jtwpa4wm','rfsquid3wm','snail3wm'})
        for name,c in data['cases'].items():
            with self.subTest(name=name):
                m=c['magnetic'];ann=c['backannotation'];l=np.array(m['inductance_H'])
                self.assertTrue(m['certification']['passed']);self.assertTrue(c['native_drc']['passed'])
                self.assertEqual(c['native_drc']['violations'],0);self.assertEqual(c['process_status'],'illustrative')
                self.assertEqual(m['source_artifact_sha256'],ann['magnetic_sha256'])
                np.testing.assert_allclose(l,np.array(ann['circuit']['inductance_H']),rtol=1e-14)
                self.assertGreater(np.linalg.eigvalsh(l)[0],0);self.assertTrue(np.any(l<0))
                scale=np.sqrt(np.outer(np.diag(l),np.diag(l)));off=~np.eye(len(l),dtype=bool);strong=off&(abs(l)/scale>=.01)
                for control,raw in c['magnetic_controls'].items():
                    self.assertEqual(m['certification']['comparison_artifacts'][control],raw['source_artifact_sha256'])
                    self.assertEqual(raw['run']['simulator']['binary_sha256'],m['run']['simulator']['binary_sha256'])
                    delta=np.array(raw['samples'][0]['L_H'])-l
                    if control!='frequency':
                        check=m['certification']['checks'][control]
                        self.assertAlmostEqual(check['max_self_relative'],max(abs(np.diag(delta)/np.diag(l))))
                        self.assertAlmostEqual(check['coupling_absolute'],max(abs(delta[off])/scale[off]))
                        self.assertLess(check['max_self_relative'],.02);self.assertLess(check['mutual_relative'],.05)
                        self.assertLess(check['coupling_absolute'],.005)
                    else:self.assertLess(m['certification']['checks'][control]['matrix_relative'],1e-7)
                for cap,e in c['electric'].items():
                    self.assertEqual(ann['electric_sha256'][cap],e['source_artifact_sha256'])
                    self.assertEqual(e['source']['gds_sha256'],ann['source']['gds_sha256'])
                    self.assertEqual(e['source']['component_id'],cap)
                    self.assertTrue(e['fit_accepted']);self.assertTrue(e['convergence']['passed'])
                    self.assertIn(1e9,e['frequency_hz']);self.assertIn(2e9,e['frequency_hz'])
                    self.assertIn(6e9,e['frequency_hz']);self.assertIn(100e9,e['frequency_hz'])
                    parent=e.get('reused_native_extraction')
                    if parent:
                        self.assertEqual(e['geometry'],parent['geometry']);self.assertEqual(e['capacitance_F'],parent['capacitance_F'])

    def test_response_against_actual_josim_and_numerical_refinements(self):
        for name,c in self.data['cases'].items():
            circuit=LayoutNetwork(c['backannotation']['circuit'])
            for label,r in c['responses'].items():
                with self.subTest(name=name,label=label):
                    d=r['drive'];p=circuit.periodic(d['pump_frequency'],d['pump_power'],flux_bias=d['flux_bias'],flux_pump=d['flux_pump'])
                    result=circuit.response(p,[d['signal_frequency']],sidebands=3);port=len(circuit.R)-1
                    s=result['waves'][0,result['bands']==0,port][0]
                    i=result['waves'][0,result['bands']==-d['mixing_order'],port][0]*np.sqrt(d['signal_frequency']/d['idler_frequency'])
                    for predicted,measured in [(abs(s)**2,r['pumped']['measured_gain']),(abs(i)**2,r['pumped']['measured_idler_conversion'])]:
                        self.assertLess(abs(10*np.log10(predicted/measured)),.25)
                    passive=circuit.passive([d['signal_frequency'],d['idler_frequency']],flux_bias=d['flux_bias'])
                    actual=np.array(r['passive']['measured_real'])+1j*np.array(r['passive']['measured_imag'])
                    self.assertLess(np.max(abs(passive-actual)),.02)
                    self.assertLess(r['pumped']['max_dc_phase_error'],1e-4);self.assertLess(p['maximum_phase'],1.)
                    for value in r['spectral_refinement'].values():self.assertLess(abs(value),.02)
                    for key in ('signal_change_db','idler_change_db'):self.assertLess(abs(r['acquisition_refinement'][key]),.02)
                    self.assertTrue(r['passed']);self.assertEqual(len(r['simulator']['binary_sha256']),64)
                    for native in c['transient_runs'][label].values():
                        self.assertEqual(hashlib.sha256(native['netlist'].encode()).hexdigest(),native['netlist_sha256'])
                        self.assertEqual(native['returncode'],0);self.assertEqual(native['analysis'],1)
                        self.assertEqual(native['corrector'],4);self.assertEqual(native['output_format'],'npy')
            self.assertTrue(c['timestep']['passed'])
            self.assertLess(max(abs(v) for v in c['timestep']['changes_db'].values()),.02)

    def test_calibrated_kit_and_rejected_attempts_remain_distinct(self):
        from twpa.em.kit_annotation import load_kit_data
        from twpa.models import TLS
        k=self.data['kit'];ann=k['backannotation'];model=load_kit_data(ann['circuit'])
        self.assertEqual(model.design.total_cells,270)
        self.assertEqual(ann['calibration']['status'],'user_supplied_not_EM_certificate')
        self.assertLess(ann['checks']['turn_sheet_relative_error'],1e-8)
        self.assertLess(k['passive_max_energy_error'],1e-10)
        r=model.coupled([k['signal_frequency_Hz']],TLS(loss_tangent=1e-3,temperature=.02),
                        frequency=np.linspace(0,8e9,801),pump_frequency=k['pump_frequency_Hz'])
        self.assertAlmostEqual(r.gain[0],k['coupled_gain'],places=9)
        self.assertGreaterEqual(r.added_noise[0],(1-1/r.gain[0])/2-1e-10)
        self.assertAlmostEqual(r.added_noise[0],k['coupled_added_noise'],places=9)
        for failed in self.data['rejected_coarse_magnetic']:
            self.assertFalse(failed['certification']['passed']);self.assertFalse(failed['certification']['checks']['contact']['passed'])
        snail=self.data['rejected_compact_snail'];self.assertFalse(snail['passed'])
        self.assertGreater(snail['checks']['coupling_absolute'],snail['coupling_gate'])
        for failed in self.data['rejected_electric'].values():
            self.assertFalse(failed['fit_accepted'])
            self.assertTrue(failed['negative_conductance_relative']>=.01 or failed['parallel_C_fit_max_relative_residual']>=.02 or failed['port_tail_relative']>=1e-3)
        weak=self.data['under_resolved_weak_idler'];self.assertFalse(weak['passed'])
        self.assertGreater(abs(weak['pumped']['idler_error_db']),.25)
        thread=self.data['dense_thread_control'];self.assertTrue(thread['passed']);self.assertLess(thread['matrix_relative'],1e-10)


if __name__=='__main__':unittest.main()
