"""Independent regression against full-sine JoSIM circuit measurements."""
import json
from pathlib import Path
import unittest
import numpy as np
from twpa.jtwpa import Chain, ResonantLoad, bloch, linear_response, four_wave_response


def chain_from(data):
    data = dict(data)
    if data['resonator'] is not None:
        data['resonator'] = ResonantLoad(**data['resonator'])
    return Chain(**data)


class JTWPAReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((Path(__file__).parent/'fixtures/josim_jtwpa.json').read_text())

    def test_passive_scattering_and_bloch(self):
        for case in self.data['cases'].values():
            d = case['summary']; c = chain_from(d['chain'])
            f = [d['signal_frequency'], d['pump_frequency'], 2*d['pump_frequency']-d['signal_frequency']]
            predicted = linear_response(c, f)
            actual = np.array(d['passive']['s21_real'])+1j*np.array(d['passive']['s21_imag'])
            self.assertLess(max(abs(actual-predicted['s21'])), .03)
            self.assertLess(max(abs(np.array(d['passive']['measured_k'])-bloch(c, f).k)), 1e-3)

    def test_four_wave_scattering_against_circuits(self):
        for name in ['bare256', 'bare384', 'rpm384']:
            d = self.data['cases'][name]['summary']; a = d['active']
            r = four_wave_response(chain_from(d['chain']), d['pump_frequency'], d['pump_power'],
                    d['signal_frequency'], d['signal_power'], harmonics=a['harmonics'],
                    sidebands=a['sidebands'], quadrature=a['quadrature'])
            self.assertLess(abs(10*np.log10(r.gain/a['measured_gain'])), .1)
            self.assertLess(abs(10*np.log10(r.idler_conversion/a['measured_conversion'])), .1)
            self.assertLess(abs(10*np.log10(r.pump_transmission/a['measured_pump_transmission'])), .1)

    def test_depletion_against_circuit(self):
        d = self.data['cases']['depleted256']['summary']; a = d['active']
        r = four_wave_response(chain_from(d['chain']), d['pump_frequency'], d['pump_power'],
                    d['signal_frequency'], d['signal_power'], depleted=True,
                    harmonics=a['harmonics'], sidebands=a['sidebands'], quadrature=a['quadrature'])
        self.assertLess(abs(10*np.log10(r.gain/a['measured_gain'])), .1)
        self.assertLess(abs(10*np.log10(r.idler_conversion/a['measured_conversion'])), .1)
        self.assertLess(abs(r.pump_transmission-a['measured_pump_transmission']), .005)
        self.assertLess(r.pump_transmission, .8)

    def test_resolution_windows_and_provenance(self):
        for name in self.data['refinement']:
            self.assertLess(max(abs(x) for x in self.data['refinement'][name].values()), .1)
        for group in self.data['model_resolution'].values():
            for row in group['results']:
                self.assertLess(abs(row['gain_change_db']), .01)
                self.assertLess(abs(row['conversion_change_db']), .01)
        windows = self.data['rpm_windows']['late_minus_early_db']
        self.assertLess(abs(windows['gain']), .03)
        self.assertLess(abs(windows['conversion']), .03)
        self.assertLess(abs(windows['pump_transmission']), .1)
        for case in self.data['cases'].values():
            self.assertEqual(len(case['hashes']), 3)
            self.assertIn('JoSIM:', case['summary']['simulator']['version'])
            self.assertIn(case['provenance']['binary_revision'], ('5f33997', '0e0ed0f'))
        self.assertFalse(self.data['cases']['bare384']['provenance']['supplied_revision_matches_binary'])
