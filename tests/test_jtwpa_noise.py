"""Quantum limits and fluctuation-dissipation, independent of JoSIM noise."""
from dataclasses import replace
import unittest
import numpy as np
from numpy.testing import assert_allclose
from scipy.constants import h, k as kb
from twpa.jtwpa import Chain, ResonantLoad, bloch, linear_response, four_wave_response, quantum_noise
from twpa.jpa_circuit import PHIBAR
from twpa.jtwpa_noise import _covariances, _quadratures


def chain(cells=16):
    return Chain(PHIBAR/2.5e-9, 1e-12, 1e-12, cells, junction_resistance=1e25)


def power(c):
    return .5*c.port_impedance*(.4*c.critical_current)**2


class FourWaveQuantumNoise(unittest.TestCase):
    def assert_canonical(self, r):
        self.assertLess(r.commutator_error, 1e-8)
        self.assertLess(r.anomalous_commutator_error, 1e-8)
        self.assertGreater(r.minimum_uncertainty_eigenvalue, -1e-8)
        assert_allclose(r.quadrature_covariance, r.quadrature_covariance.T, atol=1e-10)

    def test_ideal_amplifier_and_high_gain_limit(self):
        # Independent two-mode Bogoliubov amplifier, including G -> infinity.
        for gain in [1., 2., 10., 1e6]:
            u = np.eye(2)*np.sqrt(gain)
            v = np.array([[0., 1.], [1., 0.]])*np.sqrt(gain-1)
            normal, anomalous = _covariances(u, v, np.full(2, .5))
            added = normal[0, 0]/gain-.5
            assert_allclose(added, (1-1/gain)/2, atol=1e-15)
            assert_allclose(normal[0, 0]-.5, gain-1, rtol=1e-14)
            if gain < 100:
                # Two-mode squeezed vacuum: covariance is pure and squeezed.
                q = _quadratures(normal, anomalous)
                assert_allclose(np.linalg.det(q), 1/16, rtol=1e-10)
        self.assertAlmostEqual(added, .5, delta=1e-6)
        gain, n_idler = 5., 2.
        u = np.eye(2)*np.sqrt(gain)
        v = np.array([[0., 1.], [1., 0.]])*np.sqrt(gain-1)
        normal, _ = _covariances(u, v, np.array([.5, n_idler+.5]))
        assert_allclose(normal[0, 0]/gain-.5, (1-1/gain)*(n_idler+.5))

    def test_passive_lossless_vacuum_and_finite_scattering(self):
        c = chain()
        r = quantum_noise(c, 1e9, 0., .9e9)
        self.assert_canonical(r)
        exact = linear_response(c, r.frequency)
        for m in range(len(r.frequency)):
            assert_allclose(r.annihilation[2*m+1, 2*m], exact['s21'][m], atol=1e-12)
            assert_allclose(r.annihilation[2*m, 2*m], exact['s11'][m], atol=1e-12)
        assert_allclose(r.creation, 0, atol=1e-13)
        assert_allclose(r.covariance, np.eye(len(r.covariance))*.5, atol=1e-12)
        assert_allclose(r.added_noise, r.quantum_limit, atol=1e-12)

    def test_each_loss_family_restores_vacuum_and_commutator(self):
        base = chain()
        cases = [replace(base, junction_resistance=1000.),
                 replace(base, ground_resistance=2000.),
                 replace(base, resonator=ResonantLoad(1e-12, 1e-12, 25e-9, 1000.))]
        for c, name in zip(cases, ['junction', 'ground', 'resonator']):
            with self.subTest(bath=name):
                r = quantum_noise(c, 1e9, 0., .9e9)
                self.assert_canonical(r)
                assert_allclose(r.covariance, .5*np.eye(len(r.covariance)), atol=1e-11)
                # Attenuation without its vacuum bath violates [b,b^dagger]=1.
                incomplete = r.annihilation@r.annihilation.conj().T
                self.assertGreater(np.max(abs(incomplete-np.eye(len(incomplete)))), .01)
                probability = 1-np.sum(abs(r.annihilation[3])**2)
                assert_allclose(r.covariance_contributions[name][3, 3], probability/2, atol=1e-12)
                # Cold attenuator added noise = (1-eta)/(2*eta), including
                # reflection into the other port as an unobserved channel.
                assert_allclose(r.added_noise, (1/r.gain-1)/2, atol=1e-10)

    def test_finite_temperature_equilibrium_and_resistor_emission(self):
        c = replace(chain(), junction_resistance=1500., ground_resistance=3000.,
                    resonator=ResonantLoad(1e-12, 1e-12, 25e-9, 1200.))
        temperature = .12
        r = quantum_noise(c, 1e9, 0., .9e9, port_temperature=temperature,
                         junction_temperature=temperature, ground_temperature=temperature,
                         resonator_temperature=temperature)
        expected = np.repeat(1/np.expm1(h*r.frequency/(kb*temperature))+.5, 2)
        assert_allclose(r.covariance, np.diag(expected), atol=1e-11)
        self.assert_canonical(r)
        # Only the distributed internal circuit is hot; ports remain vacuum.
        r = quantum_noise(c, 1e9, 0., .9e9, junction_temperature=temperature,
                          ground_temperature=temperature, resonator_temperature=temperature)
        s = linear_response(c, [.9e9])
        absorption = 1-abs(s['s11'][0])**2-abs(s['s21'][0])**2
        thermal = 1/np.expm1(h*.9e9/(kb*temperature))
        assert_allclose(r.output_noise, .5+absorption*thermal, atol=1e-11)

    def test_pumped_commutators_quantum_limit_and_classical_response(self):
        c = replace(chain(64), port_impedance=float(bloch(chain(), [1e9]).impedance[0].real))
        p = power(c)
        r = quantum_noise(c, 1e9, p, .9e9)
        classical = four_wave_response(c, 1e9, p, .9e9, p*1e-7)
        self.assert_canonical(r)
        self.assertTrue(r.phase_preserving)
        assert_allclose(r.gain, classical.gain, rtol=1e-10)
        assert_allclose(abs(r.creation[5, 2])**2, classical.idler_conversion, rtol=1e-9)
        self.assertGreater(r.gain, 1.05)
        self.assertGreaterEqual(r.excess_noise, -1e-10)
        self.assertLess(r.excess_noise, .003)
        # At vacuum, the exact finite-gain excess is the unwanted
        # annihilation-channel weight, not an arbitrary noise floor.
        other_u = np.sum(abs(r.annihilation[3])**2)-r.gain+np.sum(abs(r.bath_annihilation[3])**2)
        assert_allclose(r.excess_noise, other_u/r.gain, atol=1e-11)
        self.assertLess(np.linalg.eigvalsh(r.quadrature_covariance)[0], .5)

    def test_pumped_loss_thermal_profiles_and_covariance_decomposition(self):
        c = replace(chain(32), junction_resistance=20000., ground_resistance=30000.,
                    resonator=ResonantLoad(1e-14, 1e-12, 25e-9, 100000.))
        cold = quantum_noise(c, 1e9, power(c), .9e9)
        profile = np.zeros(c.cells); profile[7] = .2
        hot = quantum_noise(c, 1e9, power(c), .9e9, junction_temperature=profile)
        self.assert_canonical(hot)
        assert_allclose(hot.gain, cold.gain, rtol=1e-13)
        self.assertGreater(hot.added_noise, cold.added_noise)
        self.assertGreaterEqual(cold.excess_noise, -1e-10)
        normal, _ = _covariances(hot.bath_annihilation, hot.bath_creation, hot.bath_occupation)
        assert_allclose(hot.covariance-cold.covariance, normal, atol=1e-12)
        assert_allclose(sum(hot.covariance_contributions.values()), hot.covariance, atol=1e-12)
        # Also exercise port temperatures by physical frequency and side.
        temp = np.zeros((len(hot.frequency), 2)); temp[2, 0] = .15
        idler_hot = quantum_noise(c, 1e9, power(c), .9e9, port_temperature=temp)
        self.assertGreater(idler_hot.added_noise, cold.added_noise)

    def test_pumped_cubic_and_resolution(self):
        c = chain(16)
        for kind in ['cubic', 'sine']:
            r = quantum_noise(c, 1e9, power(c), .9e9, nonlinearity=kind, harmonics=1)
            self.assert_canonical(r)
        a = quantum_noise(c, 1e9, power(c), .9e9)
        b = quantum_noise(c, 1e9, power(c), .9e9, harmonics=5, quadrature=24)
        assert_allclose(a.added_noise, b.added_noise, atol=2e-5)

    def test_coherently_depleted_operating_point_is_phase_sensitive(self):
        c = chain(16)
        r = quantum_noise(c, 1e9, power(c), .9e9, signal_power=power(c)*.1,
                          depleted=True, sidebands=2)
        self.assert_canonical(r)
        self.assertFalse(r.phase_preserving)
        self.assertTrue(np.isnan(r.added_noise))
        self.assertTrue(np.isnan(r.quantum_limit))
        classical = four_wave_response(c, 1e9, power(c), .9e9, power(c)*.1,
                                       depleted=True, sidebands=2)
        assert_allclose(r.node_phase, classical.node_phase, atol=1e-12)

    def test_validation_and_no_commutator_repair(self):
        c = chain(4)
        for kwargs in [dict(port_temperature=-1), dict(ground_temperature=np.nan),
                       dict(junction_temperature=[0, 1]), dict(commutator_tolerance=0)]:
            with self.assertRaises(ValueError):
                quantum_noise(c, 1e9, power(c), .9e9, **kwargs)
        with self.assertRaises(RuntimeError):
            quantum_noise(c, 1e9, power(c), .9e9, commutator_tolerance=1e-18)


if __name__ == '__main__':
    unittest.main()
