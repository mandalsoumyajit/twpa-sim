"""rf-SQUID/SNAIL potentials, 3WM conservation and quantum FDT limits."""
from dataclasses import replace
import unittest
import numpy as np
from numpy.testing import assert_allclose
from scipy.constants import h, k as kb
from twpa import three_wave as tw
from twpa.jpa_circuit import PHIBAR, Tone
from twpa.jtwpa import _basis, _jacobian, _linear_system, _real_vector, _complex_vector, _residual


def power(c, fraction=.25):
    return .5*c.port_impedance*(fraction*c.linear_chain.critical_current)**2


class ThreeWave(unittest.TestCase):
    def test_potential_equilibria_symmetry_flux_periodicity(self):
        for kind in ['rf-squid', 'snail']:
            e = tw.example_chain(kind, 4).element
            self.assertLess(abs(e.current(e.equilibrium)/e.critical_current), 1e-12)
            self.assertGreater(tw.coefficients(e)[0], 0)
            assert_allclose(tw.coefficients(e), tw.coefficients(replace(e, flux_bias=e.flux_bias+1)), atol=1e-12)
            reverse = tw.coefficients(replace(e, flux_bias=-e.flux_bias))
            assert_allclose(reverse, tw.coefficients(e)*[1, -1, 1], atol=1e-12)
            self.assertLess(abs(tw.coefficients(replace(e, flux_bias=0))[1]), 1e-14)
        with self.assertRaises(ValueError): tw.SNAIL(1e-7, .5, 3, 1e-12)
        with self.assertRaises(ValueError): tw.RFSQUID(1e-7, PHIBAR/1e-7, 1e-12)

    def test_potential_derivatives_and_effective_rc(self):
        e = tw.example_chain('snail', 4).element
        x = e.equilibrium+np.linspace(-.1, .1, 9)
        for order in [0, 1, 2]:
            numerical = (e.current(x+1e-5, order)-e.current(x-1e-5, order))/(2e-5)
            assert_allclose(numerical, e.current(x, order+1), atol=1e-16, rtol=1e-8)
        assert_allclose(e.capacitance, (.3+1/3)*e.junction_capacitance)
        assert_allclose(1/e.resistance, (.3+1/3)/e.junction_resistance)

    def test_dc_rectification_and_schur_jacobian(self):
        for kind in ['rf-squid', 'snail']:
            c = tw.example_chain(kind, 4)
            p = tw._Potential(c.element)
            idx = ((2, 0), (1, 1), (1, -1), (4, 0))
            f = [1e9, .45e9, .55e9, 2e9]
            rng = np.random.default_rng(7)
            fields = (rng.normal(size=(4, 5))+1j*rng.normal(size=(4, 5)))*.03
            branch = fields[:, :-1]-fields[:, 1:]
            waveform = np.real(_basis(idx, 24).T@branch)
            dc, current, _ = p.dc(waveform)
            self.assertGreater(np.max(abs(dc)), 1e-5)
            assert_allclose(np.mean(current, axis=0), 0, atol=3e-13)
            matrices = [_linear_system(c.linear_chain, v) for v in f]
            direction = rng.normal(size=40)
            delta = _complex_vector(direction, 5, 4); eps = 1e-7
            numerical = _real_vector((_residual(fields+eps*delta, matrices, np.zeros_like(fields), p, idx, 24)
                                      -_residual(fields-eps*delta, matrices, np.zeros_like(fields), p, idx, 24))/(2*eps))
            assert_allclose(_jacobian(fields, matrices, p, idx, 24)@direction, numerical, atol=1e-8)

    def test_passive_abcd_and_pump_off(self):
        for kind in ['rf-squid', 'snail']:
            c = tw.example_chain(kind, 8)
            r = tw.three_wave_response(c, 1e9, 0., .45e9)
            exact = tw.linear_response(c, [.45e9])
            assert_allclose(r.gain, abs(exact['s21'][0])**2, atol=1e-11)
            self.assertLess(r.idler_conversion, 1e-20)
            self.assertTrue(np.all(tw.bloch(c, [.45e9, .55e9, 1e9]).passband))

    def test_depletion_power_balance_and_three_wave_manley_rowe(self):
        for kind in ['rf-squid', 'snail']:
            base = tw.example_chain(kind, 12)
            c = replace(base, element=replace(base.element, junction_resistance=1e25))
            p = power(c, .15)
            r = tw.three_wave_response(c, 1e9, p, .45e9, p*.02, depleted=True, sidebands=3)
            assert_allclose(np.sum(r.transmitted_power+r.reflected_power), np.sum(r.incident_power), rtol=2e-8)
            # Keep only the fundamental triplet for the exact three-mode
            # photon-flux identities. All pump harmonics are projected out.
            state = tw._state(c, 1e9, p, .45e9, harmonics=1)
            _, f, idx, matrices, source, fields, _, potential = state
            source[1, 0] = 1j*np.sqrt(8*c.port_impedance*p*.02)/(c.port_impedance*c.linear_chain.critical_current)
            from twpa.jtwpa import _newton
            fields, _, _ = _newton(fields, matrices, source, nonlinearity=potential, indices=idx, quadrature=24)
            voltage = -2j*np.pi*f[:, None]*PHIBAR*fields
            inc = source[:, 0]*c.port_impedance*c.linear_chain.critical_current/2
            change = (abs(voltage[:, -1])**2+abs(voltage[:, 0]-inc)**2-abs(inc)**2)/(2*c.port_impedance*f)
            assert_allclose(change[0], -change[1], rtol=1e-5, atol=1e-31)
            assert_allclose(change[1], change[2], rtol=1e-5, atol=1e-31)

    def test_quantum_classical_response_and_limits(self):
        for kind in ['rf-squid', 'snail']:
            c = tw.example_chain(kind, 32)
            p = power(c)
            r = tw.three_wave_response(c, 1e9, p, .45e9, p*1e-7)
            q = tw.quantum_noise(c, 1e9, p, .45e9)
            assert_allclose(q.gain, r.gain, rtol=1e-10)
            assert_allclose(abs(q.creation[3, 0])**2, r.idler_conversion, rtol=1e-9)
            self.assertLess(q.commutator_error, 1e-9)
            self.assertLess(q.anomalous_commutator_error, 1e-9)
            self.assertGreater(q.minimum_uncertainty_eigenvalue, -1e-9)
            self.assertGreater(q.excess_noise, -1e-9)
            self.assertTrue(q.phase_preserving)
            self.assertEqual(q.mode_indices[0], (1, 1))

    def test_all_loss_baths_equilibrium_and_heating(self):
        from twpa.jtwpa import ResonantLoad
        for kind in ['rf-squid', 'snail']:
            c = replace(tw.example_chain(kind, 8, junction_resistance=20000., ground_resistance=1e5),
                        resonator=ResonantLoad(1e-14, 1e-12, 25e-9, 1e4))
            cold = tw.quantum_noise(c, 1e9, 0., .45e9)
            assert_allclose(cold.covariance, .5*np.eye(len(cold.covariance)), atol=1e-11)
            q = tw.quantum_noise(c, 1e9, 0., .45e9, port_temperature=.08,
                   junction_temperature=.08, ground_temperature=.08, resonator_temperature=.08)
            expected = np.repeat(1/np.expm1(h*q.frequency/(kb*.08))+.5, 2)
            assert_allclose(q.covariance, np.diag(expected), atol=1e-10)
            passive_ports = q.annihilation@q.annihilation.conj().T
            self.assertGreater(np.max(abs(passive_ports-np.eye(len(passive_ports)))), 1e-3)
            p = power(c)
            cold = tw.quantum_noise(c, 1e9, p, .45e9)
            hot = tw.quantum_noise(c, 1e9, p, .45e9, junction_temperature=.05, ground_temperature=.05)
            self.assertGreater(hot.added_noise, cold.added_noise)
            self.assertGreater(hot.excess_noise, -1e-9)
            assert_allclose(sum(hot.covariance_contributions.values()), hot.covariance, atol=1e-12)

    def test_harmonic_resolution_and_physical_netlists(self):
        for kind in ['rf-squid', 'snail']:
            c = tw.example_chain(kind, 8)
            p = power(c)
            a = tw.three_wave_response(c, 1e9, p, .45e9)
            b = tw.three_wave_response(c, 1e9, p, .45e9, harmonics=4, quadrature=32)
            assert_allclose(a.gain, b.gain, rtol=3e-4)
            deck = tw.chain_netlist(c, [Tone(1e9, p)])
            self.assertIn('KBIAS0 LLOOP0 LCOIL0', deck)
            self.assertIn('IBIAS0 0 COIL0 PWL', deck)
            self.assertEqual(deck.count('BJL0_'), 3 if kind=='snail' else 0)
            if kind=='snail':
                self.assertIn('JJS', deck)
                with self.assertRaises(ValueError): tw.chain_netlist(c, [Tone(1e9, p)], snail_loop_inductance=1e-9)
        with self.assertRaises(ValueError): tw.three_wave_response(c, 1e9, p, .5e9)


if __name__ == '__main__':
    unittest.main()
