"""Re-evaluate preserved real native matrices and circuit measurements.

These are regressions of completed field/transient evidence, not live solvers.
"""
import json
import hashlib
from pathlib import Path
import unittest
import numpy as np
from twpa.em.coupled_jpa import CoupledJPA,GalvanicJPA,pump,response


class NativeLondonReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=json.loads((Path(__file__).resolve().parents[1]/'validation/magnetic-jpa-reference.json').read_text())

    def test_signed_native_reduction_and_energy_decomposition(self):
        m=self.d['magnetic'];t=np.array(m['reduction']);l=np.array(m['inductance_H'])
        for sample in m['samples']:
            z=np.array(sample['native_Z_real_ohm'])+1j*np.array(sample['native_Z_imag_ohm'])
            np.testing.assert_allclose((t.T@z@t).imag/(2*np.pi*sample['frequency_Hz']),l,rtol=1e-12)
            self.assertLess(np.max(abs(z-z.T))/np.max(abs(z)),1e-6)
        c=self.d['magnetic_controls'];deltas=[]
        for p,n in [('plus','minus'),('plus_refined','minus_refined')]:
            delta=c[p]['settings']['lambda_scale']**2-c[n]['settings']['lambda_scale']**2
            deltas.append((np.array(c[p]['samples'][0]['L_H'])-np.array(c[n]['samples'][0]['L_H']))/delta)
        np.testing.assert_allclose(deltas[1],m['kinetic_energy_H'],rtol=1e-12)
        np.testing.assert_allclose(l-np.array(m['kinetic_energy_H']),m['magnetic_energy_H'],rtol=1e-12)
        self.assertGreater(np.linalg.eigvalsh(m['kinetic_energy_H'])[0],0)
        self.assertGreater(abs(np.array(m['SC_minus_PEC_H'])[0,0]-deltas[1][0,0]),.5e-12)
        self.assertTrue(m['certification']['passed'])
        self.assertFalse(self.d['rejected_coarse_contact']['passed'])
        self.assertGreater(self.d['rejected_coarse_contact']['mutual_relative'],.05)
        self.assertLess(m['certification']['checks']['contact']['mutual_relative'],.05)

    def test_analytic_wire_and_physical_coupling_controls(self):
        d=self.d['analytic_and_coupling_controls'];w=d['wire'];c=d['coupling']
        expected=4e-7*np.pi*(.45e-6)**2*20e-6/(2e-6*.04e-6)
        self.assertAlmostEqual(w['derivative_kinetic_H']/expected,1,places=8)
        self.assertAlmostEqual(w['lambda_squared_ratio'],4,places=6)
        sign=np.diag([1.,-1.])
        np.testing.assert_allclose(c['reversed_H'],sign@np.array(c['near_H'])@sign,rtol=1e-12)
        self.assertLess(abs(c['far_H'][0][1]),abs(c['near_H'][0][1]))

    def test_circuit_binding_harmonics_and_measured_response(self):
        r=self.d['response'];m=self.d['magnetic'];c=CoupledJPA(**r['circuit'])
        np.testing.assert_allclose(c.inductance_matrix,m['inductance_H'],rtol=1e-12)
        self.assertEqual(c.shunt_capacitance,self.d['electric']['capacitance_F'])
        self.assertAlmostEqual(c.junction_capacitance/240e-15,1)
        np.testing.assert_allclose([v['frequency'] for v in c.modes()],[v['frequency'] for v in r['modes']],rtol=1e-12)
        p=r['passive'];v=c.passive(p['measured']['signal_frequency'])
        measured=np.array(p['measured']['signal_real'])+1j*np.array(p['measured']['signal_imag'])
        self.assertLess(max(abs(v-measured)),.002)
        point=pump(c,r['pump']['frequency'],r['pump']['power'],harmonics=3)
        predicted=response(point,r['pumped']['measured']['signal_frequency'])
        self.assertLess(point.floquet_radius,1)
        for key,field in [('gain','gain'),('idler_conversion','idler_conversion')]:
            actual=np.array(r['pumped']['measured'][field])
            self.assertLess(max(abs(10*np.log10(actual/predicted[key]))),.005)
        for case in self.d['refinement_cases'].values():self.assertTrue(case['passed'])
        for runs in self.d['transient_runs'].values():
            for run in runs.values():
                self.assertEqual(hashlib.sha256(run['netlist'].encode()).hexdigest(),run['netlist_sha256'])
                self.assertEqual((run['analysis'],run['corrector'],run['output_format'],run['returncode']),(1,4,'npy',0))
        for case in self.d['sensitivity'].values():
            self.assertLess(case['signal_gain_db'],.01);self.assertLess(case['idler_conversion_db'],.01)

    def test_bare_topology_circuit_control(self):
        r=self.d['galvanic_circuit_control'];c=GalvanicJPA(**r['circuit'])
        p=pump(c,r['pump']['frequency'],r['pump']['power'],harmonics=3)
        predicted=response(p,r['pumped']['measured']['signal_frequency'])
        self.assertLess(max(abs(10*np.log10(np.array(r['pumped']['measured']['gain'])/predicted['gain']))),.005)
        self.assertTrue(r['passed'])


if __name__=='__main__':unittest.main()
