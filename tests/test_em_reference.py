"""Frozen real openEMS/JoSIM reference; optional live extraction integration."""
import json
import os
from pathlib import Path
import tempfile
import unittest
import numpy as np

from twpa import jpa
from twpa.jpa_circuit import ParallelJPA


class EMReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((Path(__file__).parents[1]/'validation/em-jpa-reference.json').read_text())

    def test_real_extraction_fit_convergence_and_scope(self):
        e=self.data['extraction']
        self.assertEqual(e['simulator']['name'],'openEMS')
        self.assertTrue(e['fit_accepted']); self.assertTrue(e['convergence']['passed'])
        self.assertEqual({c['control'] for c in e['convergence']['comparisons']},{'mesh','domain','duration'})
        self.assertLess(e['parallel_C_fit_max_relative_residual'],.02)
        self.assertLess(e['port_tail_relative'],.001)
        self.assertGreater(e['capacitance_F'],e['source']['nominal_capacitance_F']*1.1)
        self.assertIn('Interconnect',e['scope'])
        self.assertEqual(e['stack']['status'],'illustrative')

    def test_circuit_binding_and_passive_pumped_response(self):
        r=self.data['response']; c=ParallelJPA(**r['circuit']); model=c.reduced_model()
        self.assertAlmostEqual(c.shunt_capacitance,self.data['extraction']['capacitance_F'],delta=1e-27)
        self.assertAlmostEqual(c.capacitance,c.shunt_capacitance+c.junction_capacitance,delta=1e-27)
        self.assertAlmostEqual(model.frequency,r['extracted_resonance_Hz'],delta=1e-3)
        self.assertLess(r['extracted_resonance_Hz'],r['nominal_resonance_Hz'])
        self.assertLess(r['passive']['max_complex_s11_error'],.02)
        p=r['pumped']; fs=p['measured']['signal_frequency']
        point=jpa.select_pump_point(jpa.pump_points(model,p['pump_frequency_Hz'],p['pump_power_W'],phase=np.pi/2))
        predicted=jpa.kerr_response(model,point,fs)
        np.testing.assert_allclose(predicted.gain,p['model_gain'],rtol=1e-10)
        conversion=abs(predicted.scattering[:,1,0])**2
        self.assertLess(max(abs(10*np.log10(np.asarray(p['measured']['gain'])/predicted.gain))),.25)
        self.assertLess(max(abs(10*np.log10(np.asarray(p['measured']['idler_conversion'])/conversion))),.25)
        self.assertEqual(len(r['simulator_binary_sha256']),64)

    def test_refinement_exposes_failures_and_final_acceptance(self):
        cases=self.data['refinement_cases']
        self.assertTrue(cases[0]['passed'])
        self.assertTrue(cases[1]['passed'])
        self.assertFalse(cases[-1]['passed'])
        self.assertLess(self.data['sensitivity']['timestep_final_vs_fine']['signal_gain_db'],.1)
        self.assertLess(self.data['sensitivity']['timestep_final_vs_fine']['idler_conversion_db'],.1)

    @unittest.skipUnless(os.environ.get('TWPA_OPENEMS_PYTHON'),'Optional real openEMS integration')
    def test_live_gds_extraction(self):
        from twpa.em import OpenEMS,load_stack,capacitor_job
        from twpa.layout import load_process
        from twpa.layout.generate import generate,JosephsonLayout
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); gds=root/'jpa.gds'
            generate(JosephsonLayout(),load_process('jj_overlap_illustrative'),gds)
            result=OpenEMS(os.environ['TWPA_OPENEMS_PYTHON'],os.environ.get('TWPA_OPENEMS_DLL_DIRECTORY')).run(
                capacitor_job(gds,load_stack()),root/'em')
            self.assertLess(abs(result['capacitance_F']/self.data['extraction']['capacitance_F']-1),.02)


if __name__=='__main__': unittest.main()
