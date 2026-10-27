import unittest
from dataclasses import replace
import numpy as np
from numpy.testing import assert_allclose
from scipy.constants import h
from twpa.jpa import (Resonator, pump_points, select_pump_point, kerr_response,
                      flux_response, quadratures, kerr_transient)


def operating_example():
    r = Resonator(5e9, 1e7, kerr=-1e3)
    n = 3000.
    delta = 2*r.kerr*n
    fp = r.frequency+delta/(2*np.pi)
    power = n*((r.kappa/2)**2+(delta-r.kerr*n)**2)*h*fp/r.kappa_external
    point = select_pump_point(pump_points(r, fp, power))
    return r, point


class JPATests(unittest.TestCase):
    def test_passive_complex_reflection(self):
        r = Resonator(5e9, 1e7, 2e6)
        p = select_pump_point(pump_points(r, 5e9, 0.))
        f = np.array([4.995e9, 5e9, 5.004e9])
        response = kerr_response(r, p, f)
        expected = r.kappa_external/(r.kappa/2-2j*np.pi*(f-r.frequency))-1
        assert_allclose(response.scattering[:, 0, 0], expected, atol=1e-14)
        assert_allclose(response.conversion_gain, 0, atol=1e-20)
        r0 = replace(r, kappa_internal=0)
        q = kerr_response(r0, select_pump_point(pump_points(r0, 5e9, 0)), f)
        assert_allclose(q.gain, 1, atol=1e-14)

    def test_critical_coupling_and_thermal_noise(self):
        r = Resonator(5e9, 1e7, 1e7, internal_temperature=.1)
        a = flux_response(r, 10e9, 0., [5e9])
        assert_allclose(a.gain, 0, atol=1e-15)
        q = quadratures(a)
        # At critical coupling all output is the internal bath.
        occupation = 1/np.expm1(h*5e9/(1.380649e-23*.1))
        assert_allclose(q.output_covariance, np.eye(2)*(occupation+.5), atol=1e-14)

    def test_pump_cubic_and_complex_residual(self):
        r = Resonator(5e9, 1e7, 1e6, kerr=-1e3)
        fp, power, phase = 4.999e9, 8e-14, .73
        points = pump_points(r, fp, power, phase)
        self.assertGreater(len(points), 0)
        delta = 2*np.pi*(fp-r.frequency)
        incoming = np.sqrt(power/(h*fp))*np.exp(1j*phase)
        for p in points:
            assert_allclose(abs(p.amplitude)**2, p.photons, rtol=1e-10)
            residual = (1j*delta-r.kappa/2)*p.amplitude-1j*r.kerr*abs(p.amplitude)**2*p.amplitude
            assert_allclose(residual/(-np.sqrt(r.kappa_external)*incoming), 1., atol=1e-10)

    def test_bistable_roots_require_selection(self):
        r = Resonator(5e9, 1e7, kerr=-1e3)
        fp = r.frequency-r.kappa/(2*np.pi)
        power = .24*r.kappa**3*h*fp/(abs(r.kerr)*r.kappa_external)
        points = pump_points(r, fp, power)
        self.assertEqual(len(points), 3)
        self.assertEqual([p.stability for p in points], ['stable', 'unstable', 'stable'])
        with self.assertRaisesRegex(ValueError, 'multiple stable'):
            select_pump_point(points)
        low, high = select_pump_point(points, 'lower'), select_pump_point(points, 'upper')
        self.assertLess(low.photons, high.photons)
        self.assertEqual(select_pump_point(points, previous=high).photons, high.photons)
        with self.assertRaisesRegex(ValueError, 'stable'):
            kerr_response(r, points[1], [fp+.5e6])

    def test_fold_root_reported_as_marginal(self):
        r = Resonator(5e9, 1e7, kerr=-1e3)
        fp = r.frequency-r.kappa/(2*np.pi)
        delta = 2*np.pi*(fp-r.frequency)/r.kappa
        x = (-2*delta-np.sqrt(delta**2-.75))/3
        drive = x*(.25+(delta+x)**2)
        power = drive*r.kappa**3*h*fp/(abs(r.kerr)*r.kappa_external)
        points = pump_points(r, fp, power)
        marginal = [p for p in points if p.stability == 'marginal']
        self.assertGreaterEqual(len(marginal), 1)
        self.assertEqual(len(points), 2)

    def test_weak_pump_root_is_not_truncated(self):
        r = Resonator(5e9, 1e7, kerr=-1e3)
        p = pump_points(r, 5e9, 1e-35)[0]
        self.assertGreater(p.photons, 0)
        assert_allclose(p.photons, abs(p.amplitude)**2, rtol=1e-12)

    def test_zero_kerr_linear_pump(self):
        r = Resonator(5e9, 1e7)
        p = pump_points(r, 5e9, 1e-12)[0]
        assert_allclose(p.photons, 4*1e-12/(h*5e9*r.kappa_external), rtol=1e-14)
        self.assertEqual(p.stability, 'stable')

    def test_kerr_scattering_closed_form(self):
        r, p = operating_example()
        response = kerr_response(r, p, [p.frequency+.5e6, p.frequency-.5e6])
        delta = 2*np.pi*(p.frequency-r.frequency)-2*r.kerr*p.photons
        lam = -1j*r.kerr*p.amplitude**2
        omega = 2*np.pi*(response.signal_frequency-p.frequency)
        den = (r.kappa/2-1j*omega)**2+delta**2-abs(lam)**2
        expected = -1+r.kappa_external*(r.kappa/2-1j*(omega-delta))/den
        assert_allclose(response.scattering[:, 0, 0], expected, rtol=1e-13)

    def test_bogoliubov_identity_and_added_noise(self):
        r, p = operating_example()
        response = kerr_response(r, p, p.frequency+np.array([.1e6, .5e6, 1e6]))
        assert_allclose(response.gain-response.conversion_gain, 1, atol=1e-13)
        assert_allclose(response.added_noise, .5*(1-1/response.gain), atol=1e-13)

    def test_loss_ports_restore_commutator_and_heat_noise(self):
        r = Resonator(5e9, 1e7, 2e6)
        fs = [4.999e9, 5.001e9]
        cold = flux_response(r, 10e9, .3*r.kappa, fs)
        hot = flux_response(replace(r, internal_temperature=.2), 10e9, .3*r.kappa, fs)
        j = np.diag([1, -1])
        for s, b in zip(cold.scattering, cold.bath_scattering):
            assert_allclose(s@j@s.conj().T+b@j@b.conj().T, j, atol=1e-13)
        self.assertTrue(np.all(hot.added_noise > cold.added_noise))

    def test_flux_threshold_and_degenerate_quadratures(self):
        r = Resonator(5e9, 1e7)
        a = flux_response(r, 10e9, .3*r.kappa, [5e9])
        self.assertTrue(np.isnan(a.added_noise[0]))
        q = quadratures(a)
        assert_allclose(q.power_gains, [1/16, 16], atol=1e-12)
        assert_allclose(q.noise_variances, [1/32, 8], atol=1e-12)
        assert_allclose(np.linalg.det(q.output_covariance), .25, atol=1e-12)
        for lam in [.5*r.kappa, .6*r.kappa]:
            with self.assertRaisesRegex(ValueError, 'stable'):
                flux_response(r, 10e9, lam, [5e9])
        # Detuning raises the threshold, as given by |lambda|^2 < kappa^2/4+Delta^2.
        detuned = flux_response(r, 10e9+2e6, .6*r.kappa, [5.001e9])
        self.assertTrue(np.all(detuned.eigenvalues.real < 0))

    def test_phase_rotates_squeezing_axis(self):
        r = Resonator(5e9, 1e7)
        a = quadratures(flux_response(r, 10e9, .3*r.kappa, [5e9]))
        b = quadratures(flux_response(r, 10e9, .3j*r.kappa, [5e9]))
        assert_allclose(a.noise_variances, b.noise_variances, atol=1e-13)
        self.assertGreater(abs(b.output_covariance[0, 1]), 1)

    def test_transient_weak_probe_and_convergence(self):
        r, p = operating_example()
        fs = p.frequency+.5e6
        analytic = kerr_response(r, p, [fs])
        transient = kerr_transient(r, p, fs, p.power*1e-7, 20e-6, samples=2049)
        assert_allclose(transient.signal_gain, analytic.gain[0], rtol=5e-5)
        assert_allclose(transient.idler_conversion, analytic.conversion_gain[0], rtol=5e-5)
        fine = kerr_transient(r, p, fs, p.power*1e-7, 24e-6, samples=4097, rtol=1e-10, atol=1e-12)
        assert_allclose(transient.signal_gain, fine.signal_gain, rtol=5e-5)
        strong = kerr_transient(r, p, fs, p.power, 20e-6, samples=2049)
        self.assertLess(strong.signal_gain, transient.signal_gain)

    def test_input_validation(self):
        with self.assertRaises(ValueError): Resonator(5e9, -1)
        with self.assertRaises(ValueError): Resonator(5e9, 1, kerr=np.nan)
        r, p = operating_example()
        with self.assertRaises(ValueError): pump_points(r, 5e9, -1)
        with self.assertRaises(ValueError): flux_response(r, 10e9, 0, [11e9])
        with self.assertRaises(ValueError): quadratures(kerr_response(r, p, [p.frequency+1e6]))
        with self.assertRaises(ValueError): kerr_response(replace(r, frequency=6e9), p, [5e9])
        with self.assertRaises(ValueError): kerr_transient(r, p, p.frequency, 1e-15, 1e-6)


if __name__ == '__main__':
    unittest.main()
