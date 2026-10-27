"""Regression against captured outputs of the supplied original Python solvers."""
import json
from pathlib import Path
import unittest
from dataclasses import replace
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp
from twpa import preset, linear, TLS, constant_langevin, idc_fields, coupled
from twpa.nonlinear import cme_rhs

DATA = json.loads((Path(__file__).parent/'fixtures/original_solvers.json').read_text())


class OriginalRegression(unittest.TestCase):
    def test_original_coupled_gain_and_noise(self):
        d, b = preset('v2-tls')
        for row in DATA['coupled']:
            r = constant_langevin(row['kp'], row['ks'], row['ki'], row['alpha_s'],
                                 row['alpha_i'], b, d.total_cells)
            assert_allclose(r.gain, row['gain'], rtol=1e-11)
            # Original uses 2001-point trapezoidal bath integration.
            assert_allclose(r.added_noise, row['noise'], rtol=2e-6, atol=1e-9)

    def test_original_bloch_loss(self):
        d, _ = preset('v2-tls')
        d = replace(d, sections=tuple(replace(s, cell=replace(s.cell, participation=1.)) for s in d.sections))
        for row in DATA['bloch']:
            f = np.r_[0., np.array(row['frequency'])*1e9]
            # Dense grid resolves band continuation; compare exact sampled frequencies.
            grid = np.unique(np.r_[np.linspace(0, 11e9, 22001), f])
            r = linear(d, grid, tls=TLS(row['loss'], temperature=0))
            idx = np.searchsorted(grid, f[1:])
            assert_allclose(r.k[idx], row['k'], atol=1e-13)
            assert_allclose(r.alpha[idx], row['alpha'], atol=1e-13)

    def test_original_laplace(self):
        r = idc_fields()
        old = DATA['idc']
        assert_allclose(r.energy_per_length, old['energy'], rtol=2e-8)
        assert_allclose(r.participation, old['participation'], rtol=2e-8)
        assert_allclose([r.potential[i, j] for i, j in old['probes']], old['potential'], atol=1e-12)

    def test_original_nonlinear_ode(self):
        old = DATA['nonlinear']
        sol = solve_ivp(cme_rhs, (0, old['cells']), np.array(old['initial_real'], complex),
                        args=(old['epsilon'], old['k']), rtol=1e-10, atol=1e-14)
        expected = np.array(old['output_real'])+1j*np.array(old['output_imag'])
        assert_allclose(sol.y[:, -1], expected, atol=1e-13, rtol=1e-9)

    def test_ode_coupled_matches_source_transfer(self):
        d, b = preset('v2-tls')
        d = replace(d, sections=tuple(replace(s, cell=replace(s.cell, participation=1.)) for s in d.sections))
        row = DATA['coupled'][3]
        r = coupled(d, b, [4.6e9], TLS(row['loss'], temperature=0),
                    pump_frequency=DATA['pump_frequency'], pump_loss=False)
        assert_allclose(r.gain[0], row['gain'], rtol=2e-4)
        assert_allclose(r.added_noise[0], row['noise'], rtol=2e-5)


if __name__ == '__main__':
    unittest.main()
