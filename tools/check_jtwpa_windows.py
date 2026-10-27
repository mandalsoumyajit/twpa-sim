"""Check stationarity using separate late-window portions of actual JoSIM data."""
import argparse
import json
from pathlib import Path
import numpy as np
from twpa.backends.josim import load_csv, SimulationResult
from twpa.jtwpa import Chain, ResonantLoad
from twpa.josim_jtwpa import measure_chain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    summary = json.loads((args.run/'summary.json').read_text())
    parameters = dict(summary['chain'])
    if parameters['resonator'] is not None:
        parameters['resonator'] = ResonantLoad(**parameters['resonator'])
    c = Chain(**parameters)
    actual = load_csv(args.run/'probe/output.csv')
    baseline = load_csv(args.run/'pump/output.csv')
    time = next(iter(actual.values()))
    mid = (time[0]+time[-1])/2
    result = []
    fp, fs = summary['pump_frequency'], summary['signal_frequency']
    for label, select in [('early', time<=mid), ('late', time>mid)]:
        def record(data):
            return SimulationResult(args.run, {k:v[select] for k,v in data.items()}, {})
        a, b = record(actual), record(baseline)
        v = measure_chain(a, c, [fs, 2*fp-fs, fp, 3*fp], b)
        gain = float(abs(v['transmitted'][0]/v['incident'][0])**2)
        conversion = float(abs(v['transmitted'][1]/v['incident'][0])**2*fs/(2*fp-fs))
        p = measure_chain(a, c, [fp, 3*fp])
        result.append(dict(window=label, gain=gain, conversion=conversion,
                           pump_transmission=float(abs(p['transmitted'][0]/p['incident'][0])**2)))
    difference = {key:float(10*np.log10(result[1][key]/result[0][key]))
                  for key in ('gain', 'conversion', 'pump_transmission')}
    args.output.write_text(json.dumps(dict(windows=result, late_minus_early_db=difference), indent=2)+'\n')
    print(json.dumps(difference, indent=2))


if __name__ == '__main__':
    main()
