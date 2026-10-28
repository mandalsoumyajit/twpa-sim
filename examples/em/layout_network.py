"""Inspect the finite-terminal field jobs before paying for native extraction.

Run with --family squid_jpa, jtwpa4wm, rfsquid3wm, or snail3wm. This writes
checked geometry and input jobs, not accepted C/L values or fabricated gains.
See docs/layout-backannotation.md for extraction and revision decisions.
"""
import argparse
import json
from pathlib import Path
from twpa.layout import load_process
from twpa.layout.generate import JosephsonLayout,generate
from twpa.em import layout_magnetic_job,capacitor_job,load_stack


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family',choices=['squid_jpa','jtwpa4wm','rfsquid3wm','snail3wm'],default='jtwpa4wm')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cells',type=int,default=1)
    p.add_argument('--compact-snail',action='store_true',help='Two-large-junction, 40 um pitch SNAIL benchmark')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    gds=a.output/'layout.gds'
    # Geometric dimensions are micrometres. This profile is illustrative,
    # including junction densities; it is not a foundry-qualified process.
    compact=dict(large_junctions=2,cell_pitch=40,branch_height=16,capacitor_width=8,capacitor_height=8,wire_width=1,via_width=.5) if a.compact_snail else {}
    if compact and a.family!='snail3wm':raise ValueError('--compact-snail requires --family snail3wm')
    manifest=generate(JosephsonLayout(family=a.family,cells=a.cells,**compact),load_process('jj_overlap_illustrative'),gds)
    # This coarse inspection basis is intentionally inexpensive. Fixed patches
    # are physical reference planes. Refine local/bulk/thickness controls at
    # fixed patches before certifying any native FastHenry inductance matrix.
    job=layout_magnetic_job(gds,pitch_um=4,local_pitch_um=.5,nh=1,patch_side_um=.5,
                            frequencies_hz=(1e6,),max_filaments=40000)
    (a.output/'magnetic-job.json').write_text(json.dumps(job,indent=2)+'\n')
    (a.output/'input.inp').write_text(job['deck'])
    # Every CAP has its own actual coupon, including half-height end caps.
    # Translational reuse is permitted later only for exact identical coupons.
    for c in manifest['components']:
        if c['kind']!='CAP':continue
        # Include the GHz pump/probe frequencies in the electric fit band.
        electric=capacitor_job(gds,load_stack(),component_id=c['id'],
            frequencies_hz=(1e9,2e9,3e9,4e9,5e9,6e9,10e9,20e9,40e9,60e9,80e9,100e9))
        (a.output/(c['id']+'-electric-job.json')).write_text(json.dumps(electric,indent=2)+'\n')
    print('Family:',a.family,'metal terminals:',len(job['patches']),'basis branches:',len(job['ports']))
    print('Inspection filaments:',job['audit']['filaments'])
    print('Next: native extraction/refinement, strict back-annotation, then Python/JoSIM comparison.')


if __name__=='__main__':main()
