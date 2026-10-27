import unittest
from dataclasses import replace
import numpy as np
from numpy.testing import assert_allclose
from twpa.jtwpa import (Chain, ResonantLoad, bloch, linear_response, four_wave_response,
                        chain_netlist, _cube, _jacobian, _linear_system,
                        _real_vector, _complex_vector, _residual, obrien_2014, _nonlinear)
from twpa.jpa_circuit import PHIBAR, Tone, Transient


def circuit(cells=32):
    return Chain(PHIBAR/2.5e-9, 1e-12, 1e-12, cells)


class JTWPA(unittest.TestCase):
    def test_linear_bloch_analytic_and_passivity(self):
        c = circuit()
        f = np.linspace(1e6, 2e9, 301)
        w = 2*np.pi*f
        expected = np.arccos(1-w*w*c.inductance*c.ground_capacitance/(2*(1-w*w*c.inductance*c.junction_capacitance)))
        r = bloch(c, f)
        assert_allclose(r.k, expected, atol=1e-12)
        t = c.cell_matrix(f)
        assert_allclose(np.linalg.det(t), 1, atol=1e-13)
        s = linear_response(c, f)
        self.assertTrue(np.all(abs(s['s11'])**2+abs(s['s21'])**2 <= 1+1e-12))
        self.assertEqual(bloch(c, [1e9, 3e9, c.plasma_frequency]).passband.tolist(), [True, False, False])

    def test_rpm_load_and_stopband(self):
        c = replace(circuit(), resonator=ResonantLoad(10e-15, 7.036e-12, 100e-12))
        f = np.array([1e9, 5.9e9, 6e9])
        self.assertTrue(np.all(np.isfinite(linear_response(c, f)['s21'])))
        self.assertNotEqual(bloch(c, [1e9]).k[0], bloch(circuit(), [1e9]).k[0])
        paper = obrien_2014()
        self.assertEqual(paper.cells, 2000)
        self.assertAlmostEqual(paper.plasma_frequency/1e9, 27.74, delta=.1)

    def test_lossless_harmonic_power_balance(self):
        c = replace(circuit(24), junction_resistance=1e20)
        pump = .5*c.port_impedance*(.3*c.critical_current)**2
        result = four_wave_response(c, 1e9, pump, .9e9, pump*.03, depleted=True)
        assert_allclose(sum(result.transmitted_power+result.reflected_power), sum(result.incident_power), rtol=1e-8)

    def test_cube_fourier_coefficients(self):
        phasors = np.array([.3+.2j, .02-.03j, .01+.015j])[:, None]
        time = np.arange(8192)*2*np.pi/8192
        harmonics = np.array([10, 9, 11])
        waveform = (phasors[:, 0, None]*np.exp(-1j*harmonics[:, None]*time)).real.sum(axis=0)
        expected = np.array([2*np.mean(waveform**3*np.exp(1j*n*time)) for n in harmonics])
        assert_allclose(_cube(phasors)[:, 0], expected, atol=1e-14)
        sine = np.array([2*np.mean(np.sin(waveform)*np.exp(1j*n*time)) for n in harmonics])
        exact = _nonlinear(phasors, 'sine', ((1, 0), (1, 1), (1, -1)))[:, 0]
        assert_allclose(exact, phasors[:, 0]-sine, atol=1e-13)

    def test_sparse_jacobian(self):
        c = circuit(5)
        matrices = [_linear_system(c, f) for f in [1e9, .9e9, 1.1e9]]
        rng = np.random.default_rng(2)
        fields = (rng.normal(size=(3, 6))+1j*rng.normal(size=(3, 6)))*.02
        direction = rng.normal(size=36)
        eps = 1e-7
        delta = _complex_vector(direction, 6)
        numerical = _real_vector((_residual(fields+eps*delta, matrices, np.zeros_like(fields))-
                                  _residual(fields-eps*delta, matrices, np.zeros_like(fields)))/(2*eps))
        assert_allclose(_jacobian(fields, matrices)@direction, numerical, atol=1e-9, rtol=1e-8)
        indices = ((1, 0), (1, 1), (1, -1))
        numerical = _real_vector((_residual(fields+eps*delta, matrices, np.zeros_like(fields), 'sine', indices)-
                                  _residual(fields-eps*delta, matrices, np.zeros_like(fields), 'sine', indices))/(2*eps))
        assert_allclose(_jacobian(fields, matrices, 'sine', indices)@direction, numerical, atol=1e-9, rtol=1e-8)

    def test_harmonic_and_quadrature_resolution(self):
        c = circuit(32)
        pump = .5*c.port_impedance*(.3*c.critical_current)**2
        a = four_wave_response(c, 1e9, pump, .9e9, pump*1e-6, quadrature=16)
        b = four_wave_response(c, 1e9, pump, .9e9, pump*1e-6, quadrature=24)
        assert_allclose(a.gain, b.gain, rtol=1e-10)
        self.assertEqual(len(a.frequency), 6)
        self.assertIn((3, 0), a.mode_indices)
        with self.assertRaises(ValueError): four_wave_response(c, 1e9, pump, .9e9, harmonics=2)
        with self.assertRaises(ValueError): four_wave_response(c, 1e9, pump, .9e9, nonlinearity='cubic', harmonics=3)

    def test_pump_off_matches_exact_finite_chain(self):
        c = circuit()
        response = four_wave_response(c, 1e9, 0., .9e9)
        linear = linear_response(c, [.9e9])
        assert_allclose(response.gain, abs(linear['s21'][0])**2, rtol=1e-11)
        assert_allclose(response.reflection, linear['s11'][0], atol=1e-12)
        assert_allclose(response.idler_conversion, 0, atol=1e-20)

    def test_small_probe_and_manley_rowe(self):
        c = replace(circuit(64), junction_resistance=1e20)
        pump = .5*c.port_impedance*(.25*c.critical_current)**2
        a = four_wave_response(c, 1e9, pump, .9e9, pump*1e-7, harmonics=1)
        b = four_wave_response(c, 1e9, pump, .9e9, pump*1e-7, depleted=True, harmonics=1)
        assert_allclose(a.gain, b.gain, rtol=2e-5)
        strong = four_wave_response(c, 1e9, pump, .9e9, pump*.03, depleted=True, harmonics=1)
        assert_allclose(sum(strong.transmitted_power+strong.reflected_power), sum(strong.incident_power), rtol=1e-7)
        # Use frequency-normalized photon rates; h cancels from these identities.
        change = (strong.transmitted_power+strong.reflected_power-strong.incident_power)/strong.frequency
        assert_allclose(change[1], change[2], rtol=1e-6)
        assert_allclose(change[0], -2*change[2], rtol=1e-6)
        self.assertLess(strong.residual, 1e-9)

    def test_netlist_half_shunts_and_validation(self):
        c = circuit(4)
        deck = chain_netlist(c, [Tone(1e9, 1e-15)])
        self.assertIn('BJ3 N3 N4 JJMOD', deck)
        self.assertIn('CG0 N0 0 4.9999999999999999e-13', deck)
        self.assertIn('RLOAD N4 0 50', deck)
        with self.assertRaises(ValueError): replace(c, cells=0)
        with self.assertRaises(ValueError): four_wave_response(c, 1e9, 1e-15, 1e9)
        with self.assertRaises(ValueError): four_wave_response(c, 1e9, 1e-15, .3e9)


if __name__ == '__main__':
    unittest.main()
