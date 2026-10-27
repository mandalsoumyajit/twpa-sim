"""Freeze completed JoSIM measurements and numerical convergence evidence."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np


def db(values):
    return 10*np.log10(values)


def load_run(directory):
    summary = json.loads((directory/'summary.json').read_text())
    if set(summary['cases']) != {'passive', 'kerr', 'flux', 'compression'}:
        raise ValueError('capture requires the full suite including compression')
    if not all(c['passed'] for c in summary['cases'].values()):
        raise ValueError('capture refuses a failed validation suite')
    hashes = {}
    for deck in sorted(directory.glob('*/input.cir')):
        manifest = json.loads((deck.parent/'manifest.json').read_text())
        digest = hashlib.sha256(deck.read_bytes()).hexdigest()
        normalized = hashlib.sha256(deck.read_text(encoding='utf-8').encode()).hexdigest()
        if manifest['returncode'] != 0 or manifest['netlist_sha256'] not in (digest, normalized):
            raise ValueError('invalid simulator return code or netlist integrity')
        # Early Windows runs stored a hash before LF->CRLF text translation.
        # Retain both hashes; never rewrite the already-executed circuit deck.
        hashes[deck.parent.name] = dict(raw_file=digest, lf_normalized=normalized)
    return dict(summary=summary, netlist_sha256=hashes)


def differences(first, second):
    a, b = first['summary']['cases'], second['summary']['cases']
    result = {}
    for name in ('kerr', 'flux'):
        result[name+'_gain_db'] = float(max(abs(db(a[name]['measured']['gain'])-db(b[name]['measured']['gain']))))
        result[name+'_idler_db'] = float(max(abs(db(a[name]['measured']['idler_conversion'])-db(b[name]['measured']['idler_conversion']))))
    result['quadratures_db'] = float(max(abs(db(a['flux']['quadratures']['power_gains'])-db(b['flux']['quadratures']['power_gains']))))
    result['compression_gain_db'] = float(max(abs(db([r['measured_gain'] for r in a['compression']['points']])-db([r['measured_gain'] for r in b['compression']['points']]))))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True, help='Fine step, longer window')
    parser.add_argument('--fine', type=Path, required=True, help='Same fine step, standard window')
    parser.add_argument('--coarse', type=Path, required=True, help='Twice fine step, standard window')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reference, fine, coarse = map(load_run, (args.reference, args.fine, args.coarse))
    steps = [run['summary']['transient']['step'] for run in (reference, fine, coarse)]
    if steps[0] != steps[1] or not np.isclose(steps[2]/steps[1], 2):
        raise ValueError('expected equal reference/fine steps and doubled coarse step')
    payload = dict(schema_version=1, captured_utc=datetime.now(timezone.utc).isoformat(),
                   reference=reference, fine=fine, coarse=coarse,
                   timestep_change=differences(coarse, fine),
                   window_change=differences(reference, fine),
                   generator_sha256=hashlib.sha256((Path(__file__).parent/'validate_jpa_josim.py').read_bytes()).hexdigest())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2)+'\n')
    print(json.dumps(dict(timestep_change=payload['timestep_change'], window_change=payload['window_change']), indent=2))


if __name__ == '__main__':
    main()
