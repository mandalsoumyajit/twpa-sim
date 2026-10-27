import unittest
from dataclasses import replace
import numpy as np
from numpy.testing import assert_allclose
from twpa import (preset, nonlinear, small_signal, route_profile, compression,
                    idc_fields, IDCGeometry, oxide_design, coupled, TLS,
                    network_scattering, linear, Nonlinearity)


class Pipeline(unittest.TestCase):
    def test_independent_nonlinearity(self):
        d, b = preset('v2-tls')
        custom = replace(b, nonlinearity=Nonlinearity(b.epsilon, b.xi),
                         scaling_current=.003, dc_current=0)
        a = small_signal(d, b, [4.6e9])
        c = small_signal(d, custom, [4.6e9], pump_frequency=a.pump_frequency)
        assert_allclose(c.gain, a.gain, rtol=1e-12)

    def test_layout_presets_and_kerr_convention(self):
        for name, count in [('v2-tls', 49950), ('v2-folded', 41472), ('v2-four-pass', 39204)]:
            d, b = preset(name)
            self.assertEqual(d.total_cells, count)
            self.assertEqual(b.kerr_factor, 2.)
            self.assertEqual(d.sections[0].cell.capacitance, 11.16e-15)

    def test_nonlinear_small_signal_and_gradient_zero(self):
        d, b = preset('v2-tls')
        a = small_signal(d, b, [4.6e9])
        n = nonlinear(d, b, 4.6e9, b.pump_current*1e-5, pump_frequency=a.pump_frequency,
                      rtol=1e-10, atol=1e-16)
        assert_allclose(n.gain, a.gain[0], rtol=2e-5)
        z = nonlinear(d, b, 4.6e9, b.pump_current*1e-5, pump_frequency=a.pump_frequency,
                      gradient=.02, profile=route_profile('ideal'), rtol=1e-10, atol=1e-16)
        assert_allclose(z.gain, n.gain, rtol=1e-10)
        g = nonlinear(d, b, 4.6e9, b.pump_current/100, pump_frequency=a.pump_frequency,
                      gradient=.02, profile=route_profile('folded', 2))
        self.assertTrue(np.isfinite(g.gain))

    def test_pump_depletion_and_compression(self):
        d, b = preset('v2-tls')
        fp = small_signal(d, b, [4.6e9]).pump_frequency
        powers = np.linspace(-95, -55, 17)
        gains, p1db = compression(d, b, 4.6e9, powers, pump_frequency=fp)
        self.assertTrue(powers[0] < p1db < powers[-1])
        self.assertLess(gains[-1], gains[0])
        n = nonlinear(d, b, 4.6e9, 1e-6, pump_frequency=fp)
        self.assertLess(abs(n.currents[0, -1]), b.pump_current)

    def test_oxide_to_amplifier(self):
        d, b = preset('v2-tls')
        geometry = IDCGeometry(mesh=5e-9)
        coated, bare, fields = oxide_design(d, geometry, 30e-9)
        self.assertGreater(fields.energy_per_length, bare.energy_per_length)
        self.assertGreater(fields.participation, 0)
        self.assertGreater(coated.sections[0].cell.capacitance, d.sections[0].cell.capacitance)
        tls = TLS(.001, .02, 30., 2e8)
        r = coupled(coated, b, [4e9], tls)
        self.assertTrue(np.isfinite(r.gain[0]))
        self.assertGreater(r.added_noise[0], 0)
        zero, _, _ = oxide_design(d, geometry, 0.)
        self.assertEqual(zero.sections[0].cell.capacitance, d.sections[0].cell.capacitance)
        self.assertEqual(zero.sections[0].cell.participation, 0)

    def test_network_and_symmetric_cell(self):
        d, b = preset('v2-tls')
        f = np.linspace(0, 11e9, 3001)
        r = linear(d, f)
        net = network_scattering([(d.sections, d.repeats)], f)
        assert_allclose(net['s21'], r.s21, atol=1e-10)
        sym = replace(d, sections=tuple(replace(s, cell=replace(s.cell, topology='symmetric')) for s in d.sections))
        s = linear(sym, f)
        assert_allclose(s.k, r.k, atol=1e-11)
        assert_allclose(abs(s.s11)**2+abs(s.s21)**2, 1, atol=1e-8)


if __name__ == '__main__':
    unittest.main()
