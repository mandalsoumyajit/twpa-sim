"""Check harmonic, sideband and Fourier quadrature convergence for saved cases."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
from twpa.jtwpa import Chain, ResonantLoad, four_wave_response


def load_chain(parameters):
    parameters = dict(parameters)
    if parameters.get('resonator') is not None:
        parameters['resonator'] = ResonantLoad(**parameters['resonator'])
    return Chain(**parameters)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.summary.read_text())
    c = load_chain(data['chain'])
    active = data['active']
    h, b, q = active['harmonics'], active['sidebands'], active['quadrature']
    configurations = [(h, b, q), (h+2, b, q+8), (h, b, q+8)]
    if active['depleted']:
        configurations.append((h, b+1, q+8))
    results = []
    for harmonics, sidebands, quadrature in configurations:
        r = four_wave_response(c, data['pump_frequency'], data['pump_power'],
                data['signal_frequency'], data['signal_power'], depleted=active['depleted'],
                harmonics=harmonics, sidebands=sidebands, quadrature=quadrature)
        row = dict(harmonics=harmonics, sidebands=sidebands, quadrature=quadrature,
                   gain=r.gain, conversion=r.idler_conversion, pump_transmission=r.pump_transmission,
                   maximum_phase=r.maximum_junction_phase, residual=r.residual)
        results.append(row)
        print(json.dumps(row), flush=True)
    for row in results:
        row['gain_change_db'] = float(10*np.log10(row['gain']/results[0]['gain']))
        row['conversion_change_db'] = float(10*np.log10(row['conversion']/results[0]['conversion']))
        row['pump_change_db'] = float(10*np.log10(row['pump_transmission']/results[0]['pump_transmission']))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(chain=data['chain'], results=results), indent=2)+'\n')


if __name__ == '__main__':
    main()
