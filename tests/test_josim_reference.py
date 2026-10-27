"""Model regression against frozen measurements from actual JoSIM circuits."""
import json
from pathlib import Path
import unittest
import numpy as np
from numpy.testing import assert_allclose
from twpa import jpa
from twpa.jpa_circuit import ParallelJPA, SQUIDJPA


class JoSIMReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parent/'fixtures/josim_jpa.json'
        cls.data = json.loads(path.read_text())
        cls.cases = cls.data['reference']['summary']['cases']

    def assert_db(self, predicted, measured, tolerance=.25):
        error = abs(10*np.log10(predicted)-10*np.log10(measured))
        self.assertLess(float(np.max(error)), tolerance)

    def test_physical_passive_complex_reflection(self):
        row = self.cases['passive']
        r = ParallelJPA(**row['circuit']).reduced_model()
        fs = row['measured']['signal_frequency']
        response = jpa.kerr_response(r, jpa.pump_points(r, r.frequency, 0)[0], fs)
        measured = np.array(row['measured']['signal_real'])+1j*np.array(row['measured']['signal_imag'])
        self.assertLess(max(abs(response.scattering[:, 0, 0]-measured)), .02)

    def test_physical_kerr_gain_and_idler(self):
        row = self.cases['kerr']
        r = ParallelJPA(**row['circuit']).reduced_model()
        point = jpa.select_pump_point(jpa.pump_points(r, row['pump_frequency'], row['pump_power'], phase=np.pi/2))
        response = jpa.kerr_response(r, point, row['measured']['signal_frequency'])
        self.assert_db(response.gain, row['measured']['gain'])
        self.assert_db(abs(response.scattering[:, 1, 0])**2, row['measured']['idler_conversion'])

    def test_physical_kerr_compression(self):
        row = self.cases['kerr']
        timing = self.data['reference']['summary']['transient']
        r = ParallelJPA(**row['circuit']).reduced_model()
        point = jpa.select_pump_point(jpa.pump_points(r, row['pump_frequency'], row['pump_power'], phase=np.pi/2))
        for case in self.cases['compression']['points']:
            predicted = jpa.kerr_transient(r, point, row['measured']['signal_frequency'][0],
                       case['probe_power'], timing['stop'], discard=timing['start']/timing['stop'])
            self.assert_db(predicted.signal_gain, case['measured_gain'])

    def test_physical_flux_gain_quadratures_and_equilibrium(self):
        row = self.cases['flux']
        circuit = SQUIDJPA(**row['circuit'])
        r = circuit.reduced_model()
        fp, amplitude, phase = row['flux_pump']
        pair = circuit.pair_drive(amplitude, phase)
        response = jpa.flux_response(r, fp, pair, row['measured']['signal_frequency'])
        self.assert_db(response.gain, row['measured']['gain'])
        self.assert_db(abs(response.scattering[:, 1, 0])**2, row['measured']['idler_conversion'])
        q = jpa.quadratures(jpa.flux_response(r, fp, pair, [fp/2]))
        self.assert_db(q.power_gains, row['quadratures']['power_gains'])
        assert_allclose(q.external_transfer, row['quadratures']['transfer'], atol=.08)
        # A finite-window mean under the AC flux pump includes its small
        # second-order DC shift; it is not exactly the static-bias phase.
        self.assertLess(abs(circuit.junction_phase()-row['equilibrium_phase_observed']), 1e-6)

    def test_numerical_convergence_and_provenance(self):
        self.assertLess(max(self.data['timestep_change'].values()), .15)
        self.assertLess(max(self.data['window_change'].values()), .05)
        summary = self.data['reference']['summary']
        self.assertEqual(summary['analysis'], 'phase')
        self.assertIn('5f33997', summary['simulator']['version'])
        self.assertGreaterEqual(len(self.data['reference']['netlist_sha256']), 10)
