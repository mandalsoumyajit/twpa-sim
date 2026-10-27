import unittest
from dataclasses import replace
import numpy as np
from numpy.testing import assert_allclose
from scipy.constants import epsilon_0
from twpa import (Cell, Section, Design, Bias, TLS, preset, linear, cell_matrix,
                    small_signal, coupled, tls_tangent, electrostatic, default_frequency)


class Solvers(unittest.TestCase):
    def test_presets_and_lossless_energy(self):
        for name in ('v1-target', 'v1-cpw', 'v1-idc', 'v2'):
            d, b = preset(name)
            r = linear(d, default_frequency())
            assert_allclose(abs(r.s11)**2+abs(r.s21)**2, 1, atol=2e-8)
            self.assertEqual(d.total_cells, 56700 if name == 'v2' else 49361)
            self.assertEqual(r.s21[0], 1)
            g = small_signal(d, b, [3e9, 4e9])
            self.assertTrue(np.all(np.isfinite(g.gain)))

    def test_uniform_and_direct_cascade(self):
        u, l = Cell(50e-12, 10e-15), Cell(70e-12, 8e-15)
        f = np.linspace(0, 10e9, 1001)
        d = Design((Section(u, 3),), 7)
        r = linear(d, f)
        expected = np.arccos(1-u.inductance*2*u.capacitance*(2*np.pi*f)**2/2)
        assert_allclose(r.k, expected, atol=1e-12)
        d = Design((Section(u, 3), Section(l, 1)), 7, (Section(u, 2),))
        r = linear(d, f)
        for ix in (1, 500, 1000):
            a, b = cell_matrix(u, f[ix])[0], cell_matrix(l, f[ix])[0]
            t = np.linalg.matrix_power(np.linalg.matrix_power(a, 3)@b, 7)@a@a
            assert_allclose(r.s21[ix], 2/t.sum(), atol=1e-11)
            assert_allclose(r.s11[ix], (t[0, 0]+t[0, 1]-t[1, 0]-t[1, 1])/t.sum(), atol=1e-11)

    def test_grid_and_long_stopband(self):
        d, _ = preset()
        a = linear(d, default_frequency())
        b = linear(d, np.linspace(0, 11.5e9, 65537))
        assert_allclose(a.k, b.k[::2], atol=1e-11)
        c = linear(replace(d, repeats=100000), default_frequency())
        self.assertTrue(np.all(np.isfinite(c.s21_db)))

    def test_validation(self):
        with self.assertRaises(ValueError): Cell(-1, 1)
        d = Design((Section(Cell(50e-12, 10e-15), 1),), 1000)
        with self.assertRaises(ValueError): linear(d, [1., 2.])
        with self.assertRaises(ValueError): small_signal(d, Bias(), [2e9])
        r = small_signal(d, Bias(pump_current=0), [2e9, 8e9], pump_frequency=7e9)
        assert_allclose(r.gain[0], 1)
        self.assertTrue(np.isnan(r.gain[1]))

    def test_coupled_lossless_analytic_and_quantum_limit(self):
        d, b = preset()
        fs = [3e9, 4e9]
        a = small_signal(d, b, fs)
        c = coupled(d, b, fs, TLS(), pump_frequency=a.pump_frequency)
        assert_allclose(c.gain, a.gain, rtol=2e-6)
        assert_allclose(c.added_noise, (1-1/c.gain)/2, atol=2e-7)
        assert_allclose(c.pump_output_current, b.pump_current, atol=1e-14)

    def test_coupled_tls_passive_attenuator_and_noise(self):
        cell = Cell(50e-12, 10e-15, participation=.05)
        d = Design((Section(cell, 1),), 10000)
        b = Bias(pump_current=0)
        tls = TLS(loss_tangent=.01, temperature=0)
        fs, fp = 2e9, 7e9
        c = coupled(d, b, [fs], tls, pump_frequency=fp)
        r = linear(d, np.linspace(0, 10e9, 10001), tls=tls)
        alpha = np.interp(fs, r.frequency, r.alpha)
        eta = np.exp(-2*alpha*d.total_cells)
        assert_allclose(c.gain, eta, rtol=1e-6)
        assert_allclose(c.added_noise, (1-eta)/(2*eta), rtol=2e-6)
        self.assertTrue(np.all(abs(r.s11)**2+abs(r.s21)**2 <= 1+1e-9))

    def test_saturation_and_coupled_convergence(self):
        tls = TLS(.01, .02, 1., 1e8, 1e6)
        low = tls_tangent(np.array([4e9, 9e9]), tls, 0, 9e9)
        high = tls_tangent(np.array([4e9, 9e9]), tls, 1e-4, 9e9)
        self.assertLess(high[1], low[1]/1000)
        self.assertGreater(high[0], low[0]/3)
        d, b = preset()
        d = replace(d, sections=tuple(replace(s, cell=replace(s.cell, participation=.03)) for s in d.sections))
        a = small_signal(d, b, [3e9])
        c = coupled(d, b, [3e9], tls, pump_frequency=a.pump_frequency)
        e = coupled(d, b, [3e9], tls, pump_frequency=a.pump_frequency, rtol=1e-9, atol=1e-12)
        assert_allclose(c.gain, e.gain, rtol=2e-6)
        assert_allclose(c.added_noise, e.added_noise, rtol=2e-6)
        self.assertLess(c.pump_output_current[0], b.pump_current)

    def test_electrostatics_parallel_plate_and_dielectrics(self):
        er = np.ones((21, 8))*4
        fixed = np.full(er.shape, np.nan); fixed[0, :] = 0; fixed[-1, :] = 1
        r = electrostatic(er, fixed, np.ones_like(er, bool), dx=1e-8, dy=1e-8)
        assert_allclose(r.capacitance_per_length, epsilon_0*4*8/20, rtol=1e-12)
        assert_allclose(r.participation, 1, atol=1e-12)
        assert_allclose(r.potential[:, 3], np.linspace(0, 1, 21), atol=1e-12)
        # Two dielectric layers in series, checked against the face network.
        er[11:, :] = 8
        oxide = np.zeros_like(er, bool); oxide[:11, :] = True
        r = electrostatic(er, fixed, oxide, dx=1e-8, dy=1e-8, periodic_x=True)
        resistance = 10/4 + .5/4+.5/8 + 9/8
        assert_allclose(r.capacitance_per_length, epsilon_0*8/resistance, rtol=1e-12)
        assert_allclose(r.participation, (10.5/4)/resistance, rtol=1e-12)


if __name__ == '__main__':
    unittest.main()
