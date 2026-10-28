"""Compare parameterized KIT primitives with the original research generator.

Only the original get_v2_primitives definition is executed, with explicit
KLayout/NumPy/SciPy globals. Import-time file writes in the source are skipped.
Supply the original Python source explicitly; it is not a runtime dependency.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import klayout.db as db
import numpy as np
from scipy.optimize import brentq

from twpa.layout.generate import KITLayout
from twpa.layout.kit import primitives


def capture(source,output):
    tree=ast.parse(source.read_text(encoding='utf-8-sig'))
    tree.body=[f for f in tree.body if isinstance(f,ast.FunctionDef) and f.name=='get_v2_primitives']
    if len(tree.body)!=1:
        raise ValueError('Expected one original get_v2_primitives function')
    ns=dict(db=db,np=np,brentq=brentq)
    exec(compile(tree,str(source),'exec'),ns)
    ly=db.Layout(); ly.dbu=.001
    old=ns['get_v2_primitives'](ly); new=primitives(ly,KITLayout())
    result=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),primitives={},turn_max_width_um=old[-1])
    for name,a,b in zip(('straight','turn_right','turn_left','launch_left','launch_right'),old[:5],new[:5]):
        result['primitives'][name]=dict(xor_area_um2=(db.Region(a)^db.Region(b)).area()*1e-6,
                                      area_um2=db.Region(a).area()*1e-6,islands=db.Region(a).merged().count())
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if any(p['xor_area_um2'] for p in result['primitives'].values()):
        raise SystemExit('Primitive XOR differs from original')
    print('All five primitives have zero XOR difference.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--output',type=Path,default=Path('validation/layout-reference-primitives.json'))
    args=parser.parse_args(); capture(args.source,args.output)
