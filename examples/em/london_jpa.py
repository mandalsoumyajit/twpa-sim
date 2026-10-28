"""Prepare a finite-return London JPA job; no external simulator is launched.

The circuit is a resonant Josephson parametric amplifier: a junction and
shunt capacitor connected by physical metal. A neighboring pickup loop can
couple magnetic flux into that circuit. We prepare the actual conductor
current paths so an external FastHenry solve can determine their self and
signed mutual inductances, including superconducting kinetic energy.

GDS stores the drawn polygons, and its manifest records design/process
assumptions. M1/M2 are metal layers and VIA is their interlayer contact.
The stack supplies film thickness and London penetration length lambda at
a stated temperature. These illustrative material values need calibration.

The printed filaments are numerical current elements, not lithographic wires.
Mask areas are in square micrometres. Each oriented port names the two finite
electrodes defining its positive current direction. input.inp is a solver
deck and job.json records its geometry, terminal definitions and audit.
No extracted inductance, gain or quantum noise is produced by this script.

Install twpa-sim[layout]. Use tools/validate_magnetic_jpa.py to solve and
certify the field models before circuit back-annotation.
See docs/magnetic-extraction.md and the Whiteley superconducting extension:
https://wrcad.com/manual/xicmanual/node501.html.
"""
import argparse
import json
from pathlib import Path
from twpa.layout import load_process
from twpa.layout.generate import JosephsonLayout,generate
from twpa.em import jpa_magnetic_job,load_magnetic_stack


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bare',action='store_true',help='Original galvanic signal/ground reference planes')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    # Geometry units are um. This is a separate physical pickup conductor,
    # not an ideal mutual-inductance label attached to an otherwise bare mask.
    design=JosephsonLayout(inductive_coupler=not args.bare,coupler_gap=8.,coupler_height=10.,coupler_length=40.)
    gds=args.output/'jpa.gds'
    generate(design,load_process('jj_overlap_illustrative'),gds)
    # Stack lambdas are explicit operating-temperature inputs, in um.
    # pitch_um sets bulk in-plane current spacing; local_pitch_um refines
    # contacts/electrodes, while nh subdivides film thickness numerically.
    # The 1 MHz frequency is a passive London extraction sample, not the
    # nonlinear amplifier pump frequency. Start with this small inspection
    # mesh; it is NOT a convergence certificate or an inductance result.
    job=jpa_magnetic_job(gds,load_magnetic_stack(),pitch_um=4,local_pitch_um=1,nh=1,frequencies_hz=(1e6,))
    (args.output/'input.inp').write_text(job['deck'],encoding='utf-8')
    (args.output/'job.json').write_text(json.dumps(job,indent=2)+'\n',encoding='utf-8')
    print('Topology:',job['topology'],'filaments:',job['audit']['filaments'])
    print('Mask areas (um²):',job['audit']['conductor_areas_um2'])
    print('Raw ports:',[(v['name'],v['positive_patch'],v['negative_patch']) for v in job['ports']])
    print('Solve/refine this basis before using L. Total London L already includes kinetic energy.')


if __name__=='__main__':main()
