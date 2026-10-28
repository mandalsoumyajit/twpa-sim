"""Extraction contract, provenance, sizing and back-annotation regression."""
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from twpa.em import load_stack,EMStack,certify,read_extraction,back_annotate_jpa
from twpa.layout.generate import digest

try:
    import klayout.db
    HAVE_LAYOUT=True
except ImportError:
    HAVE_LAYOUT=False


def seal(d):
    d=deepcopy(d); d.pop('artifact_sha256',None); d['artifact_sha256']=digest(d)
    return d


def synthetic():
    """Synthetic fixture for contract tests; not presented as a solver result."""
    return seal(dict(schema_version=1,kind='jpa_mim_capacitor',
                     simulator=dict(name='openEMS',version='test',csxcad_version='test'),
                     source={},geometry={},stack={},stack_sha256='test',
                     settings=dict(mesh_um=.5,z_cells=3,margin_um=7,timesteps=100000),
                     run=dict(worker_sha256='test'),frequency_hz=[10e9,60e9,100e9],
                     capacitance_F=178e-15,capacitance_by_frequency_F=[178e-15]*3,
                     fit_accepted=True,convergence='unestablished'))


class EMContract(unittest.TestCase):
    def test_stack_and_rejected_inputs(self):
        s=load_stack()
        self.assertEqual(s.electrode_model,'PEC')
        for key,value in [('dielectric_thickness_um',0),('dielectric_epsilon_r',float('nan')),
                          ('electrode_model','superconducting'),('status','qualified')]:
            d=asdict(s); d[key]=value
            with self.subTest(key=key):
                with self.assertRaises(ValueError): EMStack(**d)

    def test_artifact_hash_and_fit_guard(self):
        d=synthetic(); self.assertEqual(read_extraction(d)['capacitance_F'],178e-15)
        d['capacitance_F']*=2
        with self.assertRaisesRegex(ValueError,'hash'): read_extraction(d)
        d=synthetic(); d['fit_accepted']=False
        with self.assertRaisesRegex(ValueError,'gates'): read_extraction(seal(d))

    def test_failed_runtime_records_inputs_and_preserves_existing_directory(self):
        from twpa.em import OpenEMS
        from twpa.em.openems import ExtractionError
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'run'
            runner=OpenEMS(str(Path(temp)/'no-such-python'))
            with self.assertRaises(ExtractionError): runner.run({'synthetic':'contract-test'},target)
            self.assertTrue((target/'job.json').is_file())
            self.assertTrue((target/'run.json').is_file())
            self.assertTrue((target/'stderr.log').read_text())
            with self.assertRaises(FileExistsError): runner.run({},target)

    def test_convergence_coverage_and_physical_identity(self):
        p=synthetic(); others=[]
        for key,value in [('mesh_um',1.),('margin_um',5),('timesteps',50000)]:
            r=deepcopy(p); r['settings'][key]=value; others.append(seal(r))
        certified=certify(p,others)
        self.assertTrue(certified['convergence']['passed'])
        with self.assertRaisesRegex(ValueError,'all required'): certify(p,others[:2])
        r=deepcopy(others[0]); r['stack']={'epsilon':100}
        with self.assertRaisesRegex(ValueError,'physical'): certify(p,[seal(r),*others[1:]])
        r=deepcopy(others[0]); r['settings']['margin_um']=5
        with self.assertRaisesRegex(ValueError,'only mesh'): certify(p,[seal(r),*others[1:]])
        r=deepcopy(others[0]); r['capacitance_by_frequency_F']=[200e-15]*3
        self.assertFalse(certify(p,[seal(r),*others[1:]])['convergence']['passed'])
        r=deepcopy(others[0]); r['settings']['mesh_um']=.25
        with self.assertRaisesRegex(ValueError,'finest'): certify(p,[seal(r),*others[1:]])


@unittest.skipUnless(HAVE_LAYOUT,'Install twpa-sim[layout]')
class GeometryBinding(unittest.TestCase):
    def setUp(self):
        from twpa.layout import load_process
        from twpa.layout.generate import JosephsonLayout,generate
        self.temp=tempfile.TemporaryDirectory(); self.root=Path(self.temp.name)
        self.gds=self.root/'jpa.gds'
        generate(JosephsonLayout(),load_process('jj_overlap_illustrative'),self.gds)

    def tearDown(self): self.temp.cleanup()

    def test_actual_electrodes_and_junction_area(self):
        from twpa.em import capacitor_job
        job=capacitor_job(self.gds,load_stack())
        self.assertAlmostEqual(job['geometry']['electrode_areas_um2']['M1'],196)
        # Actual upper electrode includes the clipped feed neck, not just CAP annotation.
        self.assertAlmostEqual(job['geometry']['electrode_areas_um2']['M2'],170.5)
        self.assertAlmostEqual(job['source']['nominal_capacitance_F'],144e-15)
        self.assertAlmostEqual(job['source']['junction_estimates']['Cj_F'],240e-15)
        self.assertAlmostEqual(job['source']['junction_estimates']['Ic_A'],20e-6)
        for kw in [dict(mesh_um=10),dict(z_cells=1),dict(frequencies_hz=[1,1,2])]:
            with self.assertRaises(ValueError): capacitor_job(self.gds,load_stack(),**kw)

    def test_annotation_requires_convergence_matching_gds_and_single_cj(self):
        from twpa.em import capacitor_job
        job=capacitor_job(self.gds,load_stack())
        p=synthetic(); p['source']=job['source']; p=seal(p)
        with self.assertRaisesRegex(ValueError,'certificate'): back_annotate_jpa(self.gds,p)
        others=[]
        for key,value in [('mesh_um',1.),('margin_um',5),('timesteps',50000)]:
            r=deepcopy(p); r['settings'][key]=value; others.append(seal(r))
        certified=certify(p,others)
        circuit,record=back_annotate_jpa(self.gds,certified,output=self.root/'annotated.json')
        from twpa.em import load_circuit
        self.assertEqual(circuit,load_circuit(self.root/'annotated.json'))
        self.assertAlmostEqual(circuit.shunt_capacitance,178e-15)
        self.assertAlmostEqual(circuit.capacitance,418e-15)
        self.assertAlmostEqual(circuit.junction_capacitance,240e-15)
        self.assertEqual(record['parameter_origins']['port_impedance'],'User-supplied external embedding')
        self.gds.write_bytes(self.gds.read_bytes()+b'\0')
        with self.assertRaisesRegex(ValueError,'Stale'): back_annotate_jpa(self.gds,certified)


if __name__=='__main__': unittest.main()
