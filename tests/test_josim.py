import json
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from numpy.testing import assert_allclose
from twpa.backends.josim import JoSIM, JoSIMError, load_csv, wsl_path
from twpa.jpa_circuit import ParallelJPA, SQUIDJPA, Tone, Transient, jpa_netlist, PHIBAR
from twpa.measurement import fit_tones, port_waves


class BackendTests(unittest.TestCase):
    def setUp(self):
        JoSIM.probe.cache_clear()

    def test_wsl_paths_and_arguments_preserve_spaces(self):
        self.assertEqual(wsl_path(r'C:\a folder\test'), '/mnt/c/a folder/test')
        with self.assertRaises(ValueError): wsl_path('relative/path')
        with self.assertRaises(ValueError): wsl_path(r'\\server\folder')
        cmd = JoSIM('/a path/josim-cli', 'Ubuntu')._command(r'C:\a folder\test')
        self.assertEqual(cmd[-1], '/a path/josim-cli')
        self.assertIn('/mnt/c/a folder/test', cmd)

    def test_success_manifest_and_no_overwrite(self):
        def process(command, **kwargs):
            if command[-1] == '-v':
                return subprocess.CompletedProcess(command, 0, 'JoSIM test build', '')
            if command[-1] == '-h':
                return subprocess.CompletedProcess(command, 0, '--corrector', '')
            root = kwargs['cwd']
            (root/'output.csv').write_text('time,"V(PORT)","V(SRC0)"\n0,1,2\n1,3,4\n')
            self.assertFalse(kwargs.get('shell', False))
            return subprocess.CompletedProcess(command, 0, 'completed', '')
        with tempfile.TemporaryDirectory() as tmp, patch('twpa.backends.josim.subprocess.run', side_effect=process):
            runner = JoSIM('fake path/josim', corrector=2)
            root = Path(tmp)/'new run'
            r = runner.run('* deck\n.end\n', root, metadata={'seed': 1})
            assert_allclose(r.trace('v(port)'), [1, 3])
            manifest = json.loads((root/'manifest.json').read_text())
            self.assertEqual(manifest['returncode'], 0)
            self.assertEqual(manifest['metadata']['seed'], 1)
            self.assertEqual(len(manifest['netlist_sha256']), 64)
            self.assertEqual(manifest['netlist_sha256'], hashlib.sha256((root/'input.cir').read_bytes()).hexdigest())
            self.assertIn('-c', manifest['command'])
            with self.assertRaises(FileExistsError): runner.run('new', root)

    def test_timeout_failure_missing_csv_and_capability(self):
        for failure in ('timeout', 'returncode', 'missing'):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                def process(command, **kwargs):
                    if command[-1] in ('-v', '-h'):
                        return subprocess.CompletedProcess(command, 0, 'JoSIM test build', '')
                    if failure == 'timeout':
                        raise subprocess.TimeoutExpired(command, 1, output=b'partial log', stderr=b'timeout')
                    return subprocess.CompletedProcess(command, 1 if failure == 'returncode' else 0, 'log', 'error')
                with patch('twpa.backends.josim.subprocess.run', side_effect=process):
                    root = Path(tmp)/'run'
                    with self.assertRaises(JoSIMError): JoSIM('fake', timeout=1).run('.end', root)
                    self.assertTrue((root/'manifest.json').exists())
                    self.assertTrue((root/'stderr.log').exists())
        JoSIM.probe.cache_clear()
        with patch('twpa.backends.josim.subprocess.run', return_value=subprocess.CompletedProcess([], 0, 'old upstream', '')):
            with tempfile.TemporaryDirectory() as tmp:
                with self.assertRaisesRegex(JoSIMError, 'corrector'):
                    JoSIM('upstream', corrector=4).run('.end', Path(tmp)/'run')

    def test_csv_rejects_duplicate_and_nonmonotonic(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'output.csv'
            p.write_text('time,time\n0,0\n1,1\n')
            with self.assertRaises(JoSIMError): load_csv(p)
            p.write_text('time,V\n1,0\n0,1\n')
            with self.assertRaises(JoSIMError): load_csv(p)


class CircuitMeasurementTests(unittest.TestCase):
    def test_phasor_normalization_and_port_waves(self):
        t = np.linspace(0, 2e-6, 10001)
        frequencies = np.array([5e6, 7e6])
        expected = np.array([1.2+.7j, -.4+.5j])
        y = (expected[None, :]*np.exp(-2j*np.pi*t[:, None]*frequencies)).real.sum(axis=1)+.3
        assert_allclose(fit_tones(t, y, frequencies), expected, atol=1e-13)
        ai, ao = 2.+1j, .5-.3j
        v, current = ai+ao, (ai-ao)/50
        incoming, outgoing = port_waves(v, current, 50)
        assert_allclose(incoming, ai); assert_allclose(outgoing, ao)
        with self.assertRaises(ValueError): fit_tones(t, y, [5e6, 5e6])

    def test_single_junction_mapping_and_capacitance_count(self):
        c = ParallelJPA(PHIBAR*1e-9*(2*np.pi*1e9)**2, 999e-12)
        r = c.reduced_model()
        assert_allclose(r.frequency, 1e9, rtol=1e-12)
        assert_allclose(r.kappa_external, 20e6, rtol=1e-12)
        deck = jpa_netlist(c, [Tone(1e9, 1e-18)])
        self.assertIn('BJJ PORT 0 JJMOD', deck)
        self.assertIn('RTYPE=0', deck)
        self.assertNotIn('.temp', deck.lower())
        with self.assertRaises(ValueError): jpa_netlist(c, [Tone(1e9, 0)], flux_pump=(2e9, .001, 0))

    def test_squid_equilibrium_and_bias_ramp(self):
        c = SQUIDJPA(7e-6, 999e-12)
        phase = c.junction_phase()
        assert_allclose(phase+c.arm_inductance*c.critical_current/PHIBAR*np.sin(phase), np.pi*c.flux_bias, atol=1e-12)
        opposite = SQUIDJPA(7e-6, 999e-12, flux_bias=-.2)
        assert_allclose(c.reduced_model().frequency, opposite.reduced_model().frequency, rtol=1e-12)
        self.assertLess(c.frequency_sensitivity(), 0)
        self.assertEqual(SQUIDJPA(7e-6, 999e-12, flux_bias=0).pair_drive(.001), 0)
        deck = jpa_netlist(c, [Tone(1e9, 0)], flux_pump=(2e9, .001, 0))
        self.assertIn('IBIAS 0 COIL PWL(0 0', deck)
        self.assertIn('KONE LCOIL LARM1', deck)
        self.assertIn('KTWO LCOIL LARM2 -', deck)
        with self.assertRaises(ValueError): c.pair_drive(.4)


@unittest.skipUnless(os.environ.get('TWPA_JOSIM_EXECUTABLE'), 'Set TWPA_JOSIM_EXECUTABLE for a real JoSIM integration test')
class JoSIMIntegration(unittest.TestCase):
    def test_passive_physical_resonator(self):
        from twpa import jpa
        from twpa.josim_jpa import run_circuit, sideband_measurement
        c = ParallelJPA(PHIBAR*1e-9*(2*np.pi*1e9)**2, 999e-12)
        r = c.reduced_model()
        runner = JoSIM(os.environ['TWPA_JOSIM_EXECUTABLE'], os.environ.get('TWPA_JOSIM_WSL'))
        # Integration artifacts use an explicitly selected writable parent.
        parent = os.environ.get('TWPA_JOSIM_TEST_ROOT')
        with tempfile.TemporaryDirectory(dir=parent) as tmp:
            base = run_circuit(runner, c, Path(tmp)/'base', [Tone(997e6, 0)], Transient(1e-12, 2e-6, 1e-6))
            run = run_circuit(runner, c, Path(tmp)/'probe', [Tone(1e9, 1e-22)], Transient(1e-12, 2e-6, 1e-6))
            measured = sideband_measurement(run, base, c.port_impedance, 997e6, [1e9])
            predicted = jpa.kerr_response(r, jpa.pump_points(r, 997e6, 0)[0], [1e9])
            self.assertLess(abs(measured['gain'][0]-predicted.gain[0]), .002)


if __name__ == '__main__':
    unittest.main()
