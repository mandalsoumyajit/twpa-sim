"""Generate and check every layout family; optional native KLayout Ruby DRC.

An electrical amplifier model is not yet a fabrication mask. This script
reads explicit process and geometric design profiles, constructs conductor,
junction, capacitor and contact polygons, and writes GDS layout files.
Each companion manifest records assumptions, component estimates and hashes.
The independent checker then reads the serialized geometry to test intended
connectivity and design rules (DRC: widths, spacing, enclosure and grid).

PASS means the requested geometry/checks passed for that profile; it does
not establish gain, extracted inductance or fabrication qualification.
Physical nets are connected conductor groups. Intended nets also account
for component annotations. Unassigned islands should be zero. The component
area estimates still need the field/circuit workflow for electrical accuracy.

python examples/layout/generate_all.py --output validation/runs/layout
Add --full-kit for the 49,950-cell reference, --klayout path/to/klayout for Ruby.
Outputs are regenerable GDS, manifests and geometry/DRC reports, not signoff data.
Install twpa-sim[layout]. See docs/layout.md, docs/magnetic-extraction.md and
the KLayout documentation at https://www.klayout.de/doc/manual/drc_runsets.html.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

from twpa.layout import load_process
from twpa.layout.generate import generate,load_design
from twpa.layout.native_drc import run_klayout_drc
from twpa.layout.validate import validate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('validation/runs/layout'))
    parser.add_argument('--full-kit',action='store_true')
    parser.add_argument('--klayout',help='Optional full KLayout executable for independent Ruby DRC')
    parser.add_argument('--summary',type=Path,help='Optional compact JSON evidence record')
    args=parser.parse_args()
    folder=Path(__file__).parent
    names=['kit_reference' if args.full_kit else 'kit_small','jpa','jpa_inductive','squid_jpa',
           'jtwpa4wm','rfsquid3wm','snail3wm']
    records=[]
    for name in names:
        # Choose the process explicitly. JJ values are illustrative, not vendor rules.
        process=load_process('kit_nycreates_research' if name.startswith('kit') else 'jj_overlap_illustrative')
        design=load_design(folder/(name+'.json'),process)
        path=args.output/(name+'.gds')
        manifest=generate(design,process,path)
        # Read the serialized file back through a separate geometry checker.
        result=validate(path,report=args.output/(name+'-check'))
        if not result['passed']:
            raise RuntimeError(f'{name}: {result["violations"]}')
        native=None
        if args.klayout:
            native=run_klayout_drc(path,process,args.output/(name+'-native.lyrdb'),args.klayout)
            if not native['passed']:
                raise RuntimeError(f'{name}: native DRC {native}')
        print(name, 'PASS', result['connectivity'],flush=True)
        # Component area estimates do not replace field extraction. See the
        # London JPA example for interconnect/return and pickup coupling.
        if manifest['components']:
            print('  first component:',manifest['components'][0],flush=True)
        records.append(dict(design=name,passed=result['passed'],
                            klayout_python=result['klayout_version'],
                            process_sha256=result['process_sha256'],design_sha256=result['design_sha256'],
                            gds_sha256=result['gds_sha256'],qualification=result['qualification'],
                            rule_count=len(result['rules']),connectivity=result['connectivity'],
                            components=dict(Counter(c['kind'] for c in manifest['components'])),
                            native_drc_violations=native['violations'] if native else None))
    if args.summary:
        args.summary.parent.mkdir(parents=True,exist_ok=True)
        args.summary.write_text(json.dumps(dict(schema_version=1,records=records),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
