"""Reproduce classical JPA comparisons against physical RCSJ circuits.

python tools/validate_jpa_josim.py --executable /path/to/josim-cli
    [--wsl Ubuntu] --output validation/run-name [--step 1e-12]
"""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import numpy as np
from scipy.constants import h
from twpa import jpa
from twpa.backends import JoSIM
from twpa.jpa_circuit import ParallelJPA, SQUIDJPA, Tone, Transient, PHIBAR
from twpa.josim_jpa import run_circuit, sideband_measurement, degenerate_measurement


def db(value):
    return 10*np.log10(np.asarray(value))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', default='josim-cli')
    parser.add_argument('--wsl')
    parser.add_argument('--corrector', type=int, default=0)
    parser.add_argument('--revision')
    parser.add_argument('--analysis', choices=['phase', 'voltage'], default='phase')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--step', type=float, default=1e-12)
    parser.add_argument('--stop', type=float, default=3e-6)
    parser.add_argument('--start', type=float, default=1e-6)
    parser.add_argument('--probe-power', type=float, default=1e-22)
    parser.add_argument('--compression', action='store_true', help='Also compare finite-probe Kerr transients')
    parser.add_argument('--only', choices=['all', 'passive', 'kerr', 'flux'], default='all')
    args = parser.parse_args()
    root = args.output.resolve(); root.mkdir(parents=True, exist_ok=False)
    runner = JoSIM(args.executable, args.wsl, timeout=240, corrector=args.corrector,
                   source_revision=args.revision, analysis=1 if args.analysis=='phase' else 0)
    transient = Transient(args.step, args.stop, args.start, 20e-12)
    circuit = ParallelJPA(PHIBAR*1e-9*(2*np.pi*1e9)**2, 999e-12)
    r = circuit.reduced_model()
    summary = dict(simulator=runner.probe(), analysis=args.analysis,
                   transient=asdict(transient), probe_power=args.probe_power, cases={})

    def save(label, data):
        summary['cases'][label] = data
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
        print(label, json.dumps(data), flush=True)

    if args.only in ('all', 'passive'):
        fp = r.frequency-3e6
        fs = r.frequency+np.array([-1e6, 0., 1e6])
        base = run_circuit(runner, circuit, root/'passive-base', [Tone(fp, 0)], transient)
        run = run_circuit(runner, circuit, root/'passive-probe', [Tone(f, args.probe_power) for f in fs], transient)
        measured = sideband_measurement(run, base, circuit.port_impedance, fp, fs)
        model = jpa.kerr_response(r, jpa.select_pump_point(jpa.pump_points(r, fp, 0)), fs)
        actual_s = np.array(measured['signal_real'])+1j*np.array(measured['signal_imag'])
        error = float(max(abs(actual_s-model.scattering[:, 0, 0])))
        save('passive', dict(circuit=asdict(circuit), model=asdict(r), measured=measured,
                             model_gain=model.gain.tolist(), max_complex_s11_error=error,
                             passed=error < .02))

    if args.only in ('all', 'kerr'):
        n = .28*r.kappa/abs(r.kerr)
        delta = 2*r.kerr*n
        fp = r.frequency+delta/(2*np.pi)
        power = n*((r.kappa/2)**2+(delta-r.kerr*n)**2)*h*fp/r.kappa_external
        point = jpa.select_pump_point(jpa.pump_points(r, fp, power, phase=np.pi/2))
        fs = fp+np.array([1e6, 2e6])
        tones = [Tone(fp, power)]
        base = run_circuit(runner, circuit, root/'kerr-base', tones, transient)
        run = run_circuit(runner, circuit, root/'kerr-probe', tones+[Tone(f, args.probe_power) for f in fs], transient)
        measured = sideband_measurement(run, base, circuit.port_impedance, fp, fs)
        model = jpa.kerr_response(r, point, fs)
        errors = db(measured['gain'])-model.gain_db
        conversions = abs(model.scattering[:, 1, 0])**2
        idler_errors = db(measured['idler_conversion'])-db(conversions)
        save('kerr', dict(circuit=asdict(circuit), model=asdict(r), pump_frequency=fp,
                          pump_power=power, pump_photons=point.photons, measured=measured,
                          model_gain=model.gain.tolist(), model_conversion=conversions.tolist(),
                          gain_error_db=errors.tolist(), idler_error_db=idler_errors.tolist(),
                          passed=bool(max(abs(errors))<.25 and max(abs(idler_errors))<.25)))
        if args.compression:
            rows = []
            for index, ratio in enumerate([1e-4, .03, .3]):
                ps = power*ratio
                run = run_circuit(runner, circuit, root/f'compression-{index}', tones+[Tone(fs[0], ps)], transient)
                measurement = sideband_measurement(run, base, circuit.port_impedance, fp, [fs[0]])
                predicted = jpa.kerr_transient(r, point, fs[0], ps, transient.stop,
                                               discard=transient.start/transient.stop)
                error = float(db(measurement['gain'][0])-db(predicted.signal_gain))
                rows.append(dict(probe_power=ps, measured_gain=measurement['gain'][0],
                                 model_gain=predicted.signal_gain, gain_error_db=error))
            save('compression', dict(points=rows, passed=all(abs(row['gain_error_db'])<.25 for row in rows)))

    if args.only in ('all', 'flux'):
        squid = SQUIDJPA(7e-6, 999e-12, arm_inductance=1e-12)
        sr = squid.reduced_model()
        amplitude = .5*sr.kappa/abs(squid.frequency_sensitivity())
        pair = squid.pair_drive(amplitude)
        pump = (2*sr.frequency, amplitude, 0.)
        fs = sr.frequency+np.array([1e6, 2e6])
        base = run_circuit(runner, squid, root/'flux-base', [Tone(sr.frequency, 0)], transient, flux_pump=pump)
        run = run_circuit(runner, squid, root/'flux-probe', [Tone(f, args.probe_power) for f in fs], transient, flux_pump=pump)
        measured = sideband_measurement(run, base, squid.port_impedance, sr.frequency, fs)
        model = jpa.flux_response(sr, pump[0], pair, fs)
        errors = db(measured['gain'])-model.gain_db
        conversions = abs(model.scattering[:, 1, 0])**2
        idler_errors = db(measured['idler_conversion'])-db(conversions)
        runs = [run_circuit(runner, squid, root/f'flux-quadrature-{i}',
                [Tone(sr.frequency, args.probe_power, phase)], transient, flux_pump=pump)
                for i, phase in enumerate([0., np.pi/2])]
        qm = degenerate_measurement(runs, base, squid.port_impedance, sr.frequency)
        predicted = jpa.quadratures(jpa.flux_response(sr, pump[0], pair, [sr.frequency]))
        qerrors = db(qm['power_gains'])-db(predicted.power_gains)
        transfer_error = float(np.max(abs(np.array(qm['transfer'])-predicted.external_transfer)))
        save('flux', dict(circuit=asdict(squid), model=asdict(sr), flux_pump=pump,
                          equilibrium_phase_expected=squid.junction_phase(),
                          equilibrium_phase_observed=float(np.mean(abs(base.trace('P(BJ1)')))),
                          pair_drive_real=float(pair.real), pair_drive_imag=float(pair.imag), measured=measured,
                          model_gain=model.gain.tolist(), model_conversion=conversions.tolist(),
                          gain_error_db=errors.tolist(), idler_error_db=idler_errors.tolist(),
                          quadratures=qm, model_quadrature_gains=predicted.power_gains.tolist(),
                          model_quadrature_transfer=predicted.external_transfer.tolist(),
                          max_quadrature_transfer_error=transfer_error,
                          quadrature_error_db=qerrors.tolist(),
                          passed=bool(max(abs(errors))<.25 and max(abs(idler_errors))<.25 and max(abs(qerrors))<.25 and transfer_error<.08)))
    if not all(case['passed'] for case in summary['cases'].values()):
        raise SystemExit('Circuit/model acceptance failed; inspect summary and refine/explain discrepancy.')


if __name__ == '__main__':
    main()
