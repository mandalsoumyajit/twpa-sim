"""Finite-return extraction contracts and independent nonlinear network physics."""
from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

from twpa.em.magnetic import MagneticStack,load_magnetic_stack,jpa_magnetic_job
from twpa.em.fasthenry import parse_impedance
from twpa.em.magnetic_annotation import certify_magnetic,read_magnetic
from twpa.em.coupled_jpa import CoupledJPA,GalvanicJPA,pump,response
from twpa.layout.generate import digest
from twpa.jpa_circuit import PHIBAR,Tone,Transient,jpa_netlist

try:
    import klayout.db as db
except ImportError:db=None


def seal(value):
    d=deepcopy(value);d.pop('artifact_sha256',None);d['artifact_sha256']=digest(d);return d


def synthetic(scale=1,**settings):
    """Contract-only quadratic energy fixture, not solver evidence."""
    matrix=(np.array([[95.,9.5],[9.5,48.]])+scale**2*np.array([[5.,.5],[.5,2.]]))*1e-12
    s=dict(pitch_um=1.5,local_pitch_um=.25,nh=3,patch_side_um=1.,lambda_scale=scale,frequencies_hz=[1e6]);s.update(settings)
    ports=[dict(name='signal_path',positive='nodeA',negative='nodeB',positive_patch='CAP_M2',negative_patch='JJ_M1'),
           dict(name='return_path',positive='nodeC',negative='nodeD',positive_patch='JJ_M2',negative_patch='CAP_M1'),
           dict(name='pickup',positive='nodeE',negative='nodeF',positive_patch='coil_lo',negative_patch='coil_hi')]
    return seal(dict(schema_version=1,kind='jpa_london_network',source={'test':'synthetic'},manifest_sha256='synthetic',
                     stack=asdict(MagneticStack()),settings=s,ports=ports,patches=[],reduction=[[1,0],[1,0],[0,1]],
                     samples=[dict(frequency_Hz=s['frequencies_hz'][0],L_H=matrix.tolist())],
                     run=dict(simulator={'binary_sha256':'contract-test'})))


def controls():
    return dict(bulk=synthetic(pitch_um=2),thickness=synthetic(nh=1),contact=synthetic(local_pitch_um=.5),
                frequency=synthetic(frequencies_hz=[1e11]),pec=synthetic(scale=1e-3),pec_refined=synthetic(scale=1e-4),
                plus=synthetic(scale=np.sqrt(1.1)),minus=synthetic(scale=np.sqrt(.9)),
                plus_refined=synthetic(scale=np.sqrt(1.05)),minus_refined=synthetic(scale=np.sqrt(.95)))


class MagneticContracts(unittest.TestCase):
    def test_explicit_runtime_configuration_without_legacy_path_field(self):
        from twpa.em import FastHenry
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'installation.json'
            path.write_text(json.dumps(dict(variant='Whiteley superconducting FastHenry',
                linux_binary='/opt/fasthenry',binary_sha256='test-hash',wsl_distro='Ubuntu')))
            runtime=FastHenry.from_installation(path)
            self.assertEqual(runtime.executable,'/opt/fasthenry')
            self.assertEqual(runtime.binary_sha256,'test-hash')

    def test_stack_and_physical_interfaces(self):
        s=load_magnetic_stack();z1,t1,l1=s.film('M1');z2,t2,l2=s.film('M2')
        self.assertAlmostEqual(z1+t1/2,0)
        self.assertAlmostEqual((z2-t2/2)-(z1+t1/2),.0664)
        for key,value in [('m1_lambda_um',0),('temperature_K',float('inf')),('status','qualified')]:
            d=asdict(s);d[key]=value
            with self.assertRaises(ValueError):MagneticStack(**d)

    def test_signed_matrix_parser_and_corruptions(self):
        ports=[dict(positive='N0',negative='N1'),dict(positive='N2',negative='N3')]
        text='Row 1: N0 to N1\nRow 2: N2 to N3\nImpedance matrix for frequency = 1e6 2 x 2\n0 +2j 0 -0.3j\n0 -0.3j 0 +1j\n'
        samples=parse_impedance(text,ports);self.assertEqual(samples[0][1][0,1],-.3j)
        for broken in [text.replace('N0 to N1','N1 to N0'),text.replace('0 +1j',''),text.replace('+1j','nanj'),text+'ADMITTANCE matrix']:
            with self.assertRaises(ValueError):parse_impedance(broken,ports)

    def test_energy_certificate_and_fixed_physical_terminals(self):
        p=synthetic();c=controls()
        # Refinement may renumber every native node; the finite electrode
        # identities and their orientation must stay fixed.
        c['bulk']['ports'][0]['positive']='renumbered';c['bulk']=seal(c['bulk'])
        d=certify_magnetic(p,c)
        self.assertTrue(d['certification']['passed'])
        np.testing.assert_allclose(np.array(d['kinetic_energy_H'])*1e12,[[5,.5],[.5,2]],rtol=1e-12)
        np.testing.assert_allclose(np.array(d['magnetic_energy_H'])*1e12,[[95,9.5],[9.5,48]],rtol=1e-12)
        self.assertEqual(read_magnetic(d)['artifact_sha256'],d['artifact_sha256'])
        d['inductance_H'][0][0]*=2
        with self.assertRaises(ValueError):read_magnetic(d)

    def test_reject_changed_problem_and_preserve_failed_gates(self):
        p=synthetic()
        for field in ('stack','source','patches'):
            c=controls();c['bulk'][field]={'changed':1} if field!='patches' else [dict(changed=1)];c['bulk']=seal(c['bulk'])
            with self.assertRaises(ValueError):certify_magnetic(p,c)
        c=controls();c['contact']['samples'][0]['L_H'][0][1]*=2;c['contact']=seal(c['contact'])
        self.assertFalse(certify_magnetic(p,c)['certification']['passed'])
        c=controls();c.pop('frequency')
        with self.assertRaises(ValueError):certify_magnetic(p,c)
        c=controls();c['bulk']['ports'][2]['positive_patch']='coil_hi';c['bulk']=seal(c['bulk'])
        with self.assertRaisesRegex(ValueError,'orientation'):certify_magnetic(p,c)


class CircuitPhysics(unittest.TestCase):
    def circuits(self):
        common=(20e-6,178.58e-15,240e-15)
        return [CoupledJPA(*common,135.8e-12,52.56e-12,3.588e-12),
                GalvanicJPA(*common,60e-12,50e-12,-10e-12)]

    def test_energy_balance_with_mutual_and_damping(self):
        for c in self.circuits():
            state=np.array([.13,.02,-.03,.01,.02]);v=1e-6;phi,uj,uc,i,k=state;lj=PHIBAR/c.critical_current
            gradient=np.r_[np.sin(phi),uj,c.shunt_capacitance/c.junction_capacitance*uc,
                           c.inductance_matrix@state[3:]/lj]
            derivative=float(gradient@c.rhs_scaled(state,v))
            current=c.port_incidence@state[3:];w=c.plasma_omega
            expected=(v/(PHIBAR*w)*current-c.port_impedance/(lj*w)*current**2
                      -uj**2/(c.junction_resistance*c.junction_capacitance*w)
                      -uc**2/(c.internal_resistance*c.junction_capacitance*w))
            self.assertAlmostEqual(derivative,expected,delta=1e-13)

    def test_passivity_zero_pump_and_topology(self):
        for c in self.circuits():
            f=c.modes()[0]['frequency'];fs=f*np.array([.99,1,1.01])
            self.assertLessEqual(max(abs(c.passive(fs))),1+1e-12)
            p=pump(c,f,0);r=response(p,fs[[0,2]])
            np.testing.assert_allclose(r['signal'],c.passive(fs[[0,2]]),atol=1e-12)
            np.testing.assert_allclose(r['idler_conversion'],0,atol=1e-25)
            deck=jpa_netlist(c,[Tone(f,0)],Transient(1e-14,1e-9,0,1e-13))
            self.assertIn('P(BJJ)',deck);self.assertIn('KEM',deck);self.assertEqual(deck.count('.model JJMOD'),1)

    def test_full_sine_periodic_stability_and_sideband_refinement(self):
        c=self.circuits()[0];m=c.modes()[0];fp=m['frequency']-.5*m['kappa']/(2*np.pi)
        p=pump(c,fp,5e-15,harmonics=2);q=pump(c,fp,5e-15,harmonics=3,initial=p)
        self.assertLess(p.residual,1e-8);self.assertLess(q.floquet_radius,1);self.assertLess(q.phase_peak,.6)
        fs=fp+np.array([.25,.5])*m['kappa']/(2*np.pi)
        coarse=response(p,fs,sideband_order=1)
        a=response(p,fs,sideband_order=2);b=response(q,fs,sideband_order=3)
        self.assertGreater(max(b['gain']),1.1)
        self.assertLess(max(abs(10*np.log10(a['gain']/b['gain']))),.02)
        self.assertGreater(max(abs(10*np.log10(coarse['gain']/b['gain']))),.02)

    def test_signed_mutual_and_invalid_energy(self):
        c=self.circuits()[0];d=asdict(c);d['mutual_inductance']*=-1;reverse=CoupledJPA(**d)
        fs=np.linspace(c.modes()[0]['frequency']*.99,c.modes()[0]['frequency']*1.01,5)
        np.testing.assert_allclose(c.passive(fs),reverse.passive(fs),atol=1e-12)
        self.assertIn(' -',reverse.components()[-1])
        d['mutual_inductance']=1e-9
        with self.assertRaises(ValueError):CoupledJPA(**d)

    def test_reject_nonfinite_and_aliased_response_inputs(self):
        c=self.circuits()[0];f=c.modes()[0]['frequency'];p=pump(c,f,0)
        for fs in ([float('nan')],[],[[f]],[-f],[2*f]):
            with self.assertRaises(ValueError):response(p,fs)
        for kwargs in (dict(quadrature=16),dict(quadrature=128.5),dict(sideband_order=True)):
            with self.assertRaises(ValueError):response(p,[f*.99],**kwargs)
        with self.assertRaises(ValueError):pump(c,f,0,harmonics=100,quadrature=64)


@unittest.skipIf(db is None,'Install twpa-sim[layout]')
class ActualMaskBasis(unittest.TestCase):
    def test_network_annotation_binding_and_no_double_added_kinetic_energy(self):
        from twpa.layout import load_process
        from twpa.layout.generate import JosephsonLayout,generate
        from twpa.em import capacitor_job,load_stack,back_annotate_network_jpa,load_circuit
        from twpa.em.backannotation import certify
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for pickup in (False,True):
                gds=root/('pickup.gds' if pickup else 'bare.gds')
                generate(JosephsonLayout(inductive_coupler=pickup),load_process('jj_overlap_illustrative'),gds)
                job=jpa_magnetic_job(gds,pitch_um=4,local_pitch_um=1,nh=1,frequencies_hz=(1e6,))
                # Synthetic contracts bound to genuine geometry: no native field
                # accuracy is claimed by this annotation/staleness test.
                magnetic=certify_magnetic(synthetic(),controls())
                magnetic.update(source=job['source'],manifest_sha256=job['manifest_sha256'],topology=job['topology'])
                magnetic=seal(magnetic)
                from test_em import synthetic as electric_fixture
                electric=electric_fixture();ejob=capacitor_job(gds,load_stack())
                electric.update(source=ejob['source'],stack=ejob['stack']);electric=seal(electric)
                refinements=[]
                for key,value in [('mesh_um',1),('margin_um',5),('timesteps',50000)]:
                    d=deepcopy(electric);d['settings'][key]=value;refinements.append(seal(d))
                electric=certify(electric,refinements)
                circuit,record=back_annotate_network_jpa(gds,electric,magnetic,output=root/'circuit.json')
                self.assertIsInstance(circuit,CoupledJPA if pickup else GalvanicJPA)
                np.testing.assert_allclose(circuit.inductance_matrix,magnetic['inductance_H'],rtol=1e-12)
                self.assertEqual(circuit,load_circuit(root/'circuit.json'))
                self.assertAlmostEqual(circuit.junction_capacitance/240e-15,1)
                d=deepcopy(magnetic);d['certification']['passed']=False
                with self.assertRaisesRegex(ValueError,'Certified'):back_annotate_network_jpa(gds,electric,seal(d))
                d=deepcopy(magnetic);d['stack']['dielectric_gap_um']*=2
                with self.assertRaisesRegex(ValueError,'gap'):back_annotate_network_jpa(gds,electric,seal(d))
                manifest=gds.with_suffix('.manifest.json');manifest.write_bytes(manifest.read_bytes()+b' ')
                with self.assertRaisesRegex(ValueError,'manifest'):back_annotate_network_jpa(gds,electric,magnetic)
                # Restore bytes for a separate actual-GDS staleness check.
                manifest.write_bytes(manifest.read_bytes()[:-1]);gds.write_bytes(gds.read_bytes()+b'\0')
                with self.assertRaisesRegex(ValueError,'Stale'):back_annotate_network_jpa(gds,electric,magnetic)

    def test_full_conductor_coverage_ports_and_vias_for_both_topologies(self):
        from twpa.layout import load_process
        from twpa.layout.generate import JosephsonLayout,generate
        from twpa.layout.validate import validate
        with tempfile.TemporaryDirectory() as temp:
            for pickup in (False,True):
                gds=Path(temp)/('pickup.gds' if pickup else 'bare.gds')
                generate(JosephsonLayout(inductive_coupler=pickup),load_process('jj_overlap_illustrative'),gds)
                self.assertTrue(validate(gds)['passed'])
                job=jpa_magnetic_job(gds,pitch_um=4,local_pitch_um=1,nh=1,frequencies_hz=(1e6,))
                self.assertEqual(job['topology'],'pickup' if pickup else 'galvanic')
                self.assertEqual(len(job['ports']),3 if pickup else 4)
                self.assertEqual(job['audit']['conductor_components'],3 if pickup else 2)
                self.assertAlmostEqual(job['audit']['conductor_areas_um2']['M1'],814 if pickup else 622)
                self.assertEqual(len(job['audit']['vias_um']),2)
                pairs=[tuple(sorted((s['a'],s['b']))) for s in job['segments']]
                self.assertEqual(len(pairs),len(set(pairs)))
                self.assertTrue(any(s['role']=='via_film_or_fill' for s in job['segments']))
                for bad in (dict(lambda_scale=0),dict(patch_side_um=3),dict(nh=0)):
                    with self.assertRaises(ValueError):jpa_magnetic_job(gds,**bad)


if __name__=='__main__':unittest.main()
