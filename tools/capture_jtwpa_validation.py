"""Freeze actual Josephson-chain measurements and convergence evidence."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import re
from pathlib import Path
import numpy as np


def load_run(path):
    summary = json.loads((path/'summary.json').read_text())
    if not summary['passed']:
        raise ValueError(f'Cannot capture failed validation: {path}')
    hashes = {}
    for deck in sorted(path.glob('*/input.cir')):
        manifest = json.loads((deck.parent/'manifest.json').read_text())
        digest = hashlib.sha256(deck.read_bytes()).hexdigest()
        if manifest['returncode'] != 0 or digest != manifest['netlist_sha256']:
            raise ValueError('Simulator failure or netlist integrity mismatch')
        hashes[deck.parent.name] = dict(netlist_sha256=digest,
                output_sha256=hashlib.sha256((deck.parent/'output.csv').read_bytes()).hexdigest())
    simulator = summary['simulator']
    revision = re.search(r'v2\.7\.([0-9a-f]{7})', simulator['version'])
    reported = revision.group(1) if revision else None
    return dict(summary=summary, hashes=hashes,
                provenance=dict(binary_revision=reported,
                    supplied_revision_matches_binary=bool(reported and simulator.get('source_revision', '').startswith(reported))))


def changes(a, b):
    return {key:float(10*np.log10(a['summary']['active'][field]/b['summary']['active'][field]))
            for key, field in [('gain_db', 'measured_gain'), ('conversion_db', 'measured_conversion'),
                               ('pump_db', 'measured_pump_transmission')]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    names = dict(bare256='bare256fine', bare384='bare384fine',
                 depleted256='bare256depletedfinal', rpm384='rpm384fine')
    cases = {name:load_run(args.root/path) for name, path in names.items()}
    coarse = dict(bare256=load_run(args.root/'bare256'), rpm384=load_run(args.root/'rpm384steady'))
    resolutions = {name:json.loads((args.root/filename).read_text()) for name, filename in
                   [('bare256','bare-resolution.json'), ('rpm384','rpm-resolution.json'),
                    ('depleted256','depleted-resolution.json')]}
    windows = json.loads((args.root/'rpm-windows.json').read_text())
    payload = dict(schema_version=1, captured_utc=datetime.now(timezone.utc).isoformat(),
                    cases=cases, coarse=coarse, model_resolution=resolutions,
                    refinement={name:changes(cases[name], coarse[name]) for name in coarse},
                    rpm_windows=windows,
                    generator_sha256=hashlib.sha256((Path(__file__).parent/'validate_jtwpa_josim.py').read_bytes()).hexdigest())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2)+'\n')
    print(json.dumps(payload['refinement'], indent=2))


if __name__ == '__main__':
    main()
