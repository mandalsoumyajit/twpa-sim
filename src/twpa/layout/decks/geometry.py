"""KLayout Python geometry deck (standalone klayout.db, no GUI required).

python path/to/geometry.py input.gds --manifest input.manifest.json --report check
Checks exported polygons, circuit annotations, conductor connectivity and DRC.
Use the independent process.drc runset for the Ruby DRC engine as well.
"""
import argparse
import json
from pathlib import Path

from twpa.layout.validate import validate

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('input',type=Path)
parser.add_argument('--manifest',type=Path)
parser.add_argument('--report',type=Path,required=True)
args=parser.parse_args()
result=validate(args.input,args.manifest,args.report)
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] else 1)
