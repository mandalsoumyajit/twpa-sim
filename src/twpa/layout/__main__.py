"""Generate or check a mask: python -m twpa.layout --help."""
import argparse
import json
from pathlib import Path

from .generate import generate, load_design
from .process import load_process
from .validate import validate


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    g=sub.add_parser('generate',help='Write GDS and process/design manifest')
    g.add_argument('--process',required=True,help='Bundled name or local JSON path')
    g.add_argument('--design',type=Path,required=True)
    g.add_argument('--output',type=Path,required=True)
    v=sub.add_parser('check',help='Check exported GDS, geometry, connectivity and process DRC')
    v.add_argument('gds',type=Path)
    v.add_argument('--manifest',type=Path)
    v.add_argument('--report',type=Path,required=True,help='Output stem (.json and .lyrdb)')
    v.add_argument('--geometry-only',action='store_true')
    args=p.parse_args()
    if args.command=='generate':
        process=load_process(args.process)
        m=generate(load_design(args.design,process),process,args.output)
        print(json.dumps(dict(gds=str(args.output),family=m['family'],process=m['process']['id'])))
    else:
        r=validate(args.gds,args.manifest,args.report,drc=not args.geometry_only)
        print(json.dumps(dict(passed=r['passed'],violations=r['violations'],connectivity=r['connectivity'])))
        return 0 if r['passed'] else 1
    return 0


if __name__=='__main__':
    raise SystemExit(main())
