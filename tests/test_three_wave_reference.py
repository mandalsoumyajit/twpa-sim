"""Regression against physical JoSIM flux-biased junction-loop circuits."""
import json
from pathlib import Path
import unittest
import numpy as np
from twpa import three_wave as tw
from twpa.jtwpa import ResonantLoad


def chain_from(summary):
    data=dict(summary['chain'])
    cls=tw.RFSQUID if summary['kind']=='rf-squid' else tw.SNAIL
    data['element']=cls(**data['element'])
    if data['resonator'] is not None:
        data['resonator']=ResonantLoad(**data['resonator'])
    return tw.ThreeWaveChain(**data)


class ThreeWaveReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((Path(__file__).parent/'fixtures/josim_three_wave.json').read_text())

    def test_passive_physical_dispersion_and_scattering(self):
        for case in self.data['cases'].values():
            d=case['summary'];c=chain_from(d)
            f=[d['signal_frequency'],d['pump_frequency'],d['pump_frequency']-d['signal_frequency']]
            s=tw.linear_response(c,f)
            measured=np.array(d['passive']['s21_real'])+1j*np.array(d['passive']['s21_imag'])
            self.assertLess(max(abs(measured-s['s21'])),.03)
            self.assertLess(max(abs(np.array(d['passive']['measured_k'])-tw.bloch(c,f).k)),1e-3)

    def test_gain_conversion_dc_and_loss_depletion(self):
        for name,case in self.data['cases'].items():
            d=case['summary'];c=chain_from(d);a=d['active']
            r=tw.three_wave_response(c,d['pump_frequency'],d['pump_power'],d['signal_frequency'],
                         d['signal_power'],depleted=a['depleted'],**d['controls'])
            self.assertLess(abs(10*np.log10(r.gain/a['measured_gain'])),.1)
            self.assertLess(abs(10*np.log10(r.idler_conversion/a['measured_conversion'])),.1)
            self.assertLess(abs(10*np.log10(r.pump_transmission/a['measured_pump_transmission'])),.1)
            model_dc=c.element.equilibrium+r.dc_phase[max(1,c.cells//2)]
            self.assertLess(abs(model_dc-a['measured_mean_phase']),.003)
            if a['depleted']:
                self.assertLess(r.pump_transmission,.5)

    def test_convergence_geometry_and_execution_provenance(self):
        for case in self.data['cases'].values():
            artifacts=case['artifacts']
            self.assertEqual(artifacts['pump']['execution_version'],artifacts['probe']['execution_version'])
            for row in artifacts.values():
                self.assertEqual(len(row['netlist_sha256']),64)
                self.assertEqual(len(row['output_sha256']),64)
        rows=self.data['geometry']['rows']
        self.assertLess(abs(rows[1]['equilibrium_error']),abs(rows[0]['equilibrium_error'])*.6)
        for row in self.data['convergence']['depleted_resolution']:
            for key in ['gain_change_db','conversion_change_db','pump_change_db']:
                self.assertLess(abs(row[key]),.01)
        for name,data in self.data['convergence']['windows'].items():
            self.assertLess(abs(data['late_minus_early_db']['gain']),.03)
            self.assertLess(abs(data['late_minus_early_db']['conversion']),.04)
            self.assertLess(abs(data['late_minus_early_db']['pump']),.2)


if __name__=='__main__':
    unittest.main()
