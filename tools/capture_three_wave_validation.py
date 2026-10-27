"""Freeze actual 3WM circuit measurements, per-execution builds and run hashes."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def artifacts(directory):
    manifest=json.loads((directory/'manifest.json').read_text())
    deck=(directory/'input.cir').read_bytes()
    digest=hashlib.sha256(deck).hexdigest()
    if manifest['returncode']!=0 or digest!=manifest['netlist_sha256']:
        raise ValueError(f'Failed simulation or altered netlist: {directory}')
    banner=[line for line in (directory/'stdout.log').read_text().splitlines() if line.startswith('v2.')]
    if len(banner)!=1:
        raise ValueError('Missing actual execution version banner')
    return dict(netlist_sha256=digest,output_sha256=hashlib.sha256((directory/'output.csv').read_bytes()).hexdigest(),
        execution_version=banner[0],cached_probe_version=manifest['simulator']['version'].splitlines()[-1],
        physical_junctions=sum(line.startswith(b'B') for line in deck.splitlines()),
        directory=directory.name)


def load_run(root):
    summary=json.loads((root/'summary.json').read_text())
    if not summary['passed']:
        raise ValueError(f'Cannot capture failed run {root}')
    hashes={name:artifacts(root/(summary.get('baseline_directory','pump') if name=='pump' else name))
            for name in ['passive','pump','probe']}
    if hashes['pump']['execution_version']!=hashes['probe']['execution_version']:
        raise ValueError('Pump and probe actual execution builds differ')
    count=summary['chain']['cells']*(1 if summary['kind']=='rf-squid' else summary['chain']['element']['junctions']+1)
    if any(row['physical_junctions']!=count for row in hashes.values()):
        raise ValueError('Circuit did not contain the expected physical junctions')
    return dict(summary=summary,artifacts=hashes)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    cases={name:load_run(args.root/name) for name in
           ['rf128-fine','snail128-fine','rf128-loss','snail128-loss-depleted']}
    coarse={name:load_run(args.root/name) for name in ['rf128-smooth','snail128-smooth']}
    geometry=json.loads((args.root/'geometry/summary.json').read_text())
    geometry['artifacts']={name:artifacts(args.root/'geometry'/name) for name in ['nominal','half']}
    payload=dict(kind='Classical JoSIM measurements; quantum fields in summaries are model calculations',
        captured_utc=datetime.now(timezone.utc).isoformat(),cases=cases,coarse=coarse,geometry=geometry,
        convergence=json.loads((args.root/'convergence.json').read_text()))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2)+'\n')
    print('Captured',len(cases),'physical 3WM benchmarks')


if __name__=='__main__':
    main()
