"""Physical 4WM JTWPA comparisons, preserving JoSIM decks and traces."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import numpy as np
from twpa.backends import JoSIM
from twpa.jtwpa import Chain, ResonantLoad, bloch, linear_response, four_wave_response
from twpa.jpa_circuit import PHIBAR, Tone, Transient
from twpa.josim_jtwpa import run_chain, measure_chain, measure_bloch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', default='josim-cli')
    parser.add_argument('--wsl')
    parser.add_argument('--corrector', type=int, default=0)
    parser.add_argument('--revision')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cells', type=int, default=128)
    parser.add_argument('--step', type=float, default=2e-12)
    parser.add_argument('--stop', type=float, default=200e-9)
    parser.add_argument('--start', type=float, default=100e-9)
    parser.add_argument('--pump-ratio', type=float, default=.4)
    parser.add_argument('--probe-ratio', type=float, default=1e-6)
    parser.add_argument('--signal', type=float, default=.9e9)
    parser.add_argument('--rpm-frequency', type=float)
    parser.add_argument('--resonator-capacitance', type=float, default=1e-12)
    parser.add_argument('--coupling-capacitance', type=float, default=.01e-12)
    parser.add_argument('--resonator-resistance', type=float, default=1e9)
    parser.add_argument('--depleted', action='store_true')
    parser.add_argument('--harmonics', type=int, default=3)
    parser.add_argument('--sidebands', type=int, default=1)
    parser.add_argument('--quadrature', type=int, default=16)
    parser.add_argument('--ramp', type=float, default=0., help='Source amplitude turn-on time constant')
    args = parser.parse_args()
    root = args.output.resolve(); root.mkdir(parents=True, exist_ok=False)
    runner = JoSIM(args.executable, args.wsl, timeout=360, corrector=args.corrector,
                   source_revision=args.revision)
    chain = Chain(PHIBAR/2.5e-9, 1e-12, 1e-12, args.cells)
    if args.rpm_frequency is not None:
        cr = args.resonator_capacitance
        chain = replace(chain, resonator=ResonantLoad(args.coupling_capacitance, cr,
                        1/((2*np.pi*args.rpm_frequency)**2*cr), args.resonator_resistance))
    fp = 1e9
    chain = replace(chain, port_impedance=float(bloch(chain, [fp]).impedance[0].real))
    transient = Transient(args.step, args.stop, args.start, 10e-12)
    power = .5*chain.port_impedance*(args.pump_ratio*chain.critical_current)**2
    probe = power*args.probe_ratio
    fs = args.signal; fi = 2*fp-fs
    result = dict(chain=asdict(chain), transient=asdict(transient), simulator=runner.probe(),
                  pump_frequency=fp, pump_power=power, signal_frequency=fs, signal_power=probe,
                  source_ramp=args.ramp)
    def save():
        (root/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    frequencies = np.array([fs, fp, fi])
    passive_power = min(probe, power*1e-6)
    result['passive_probe_power'] = passive_power
    passive = run_chain(runner, chain, root/'passive', [Tone(f, passive_power, ramp_time=args.ramp) for f in frequencies], transient)
    measured = measure_chain(passive, chain, frequencies)
    exact = linear_response(chain, frequencies)
    s21 = measured['transmitted']/measured['incident']
    s11 = measured['reflected']/measured['incident']
    theta = measure_bloch(passive, chain, frequencies)
    predicted_bloch = bloch(chain, frequencies)
    result['passive'] = dict(s21_real=s21.real.tolist(), s21_imag=s21.imag.tolist(),
                             s11_real=s11.real.tolist(), s11_imag=s11.imag.tolist(),
                             max_s21_error=float(max(abs(s21-exact['s21']))),
                             max_s11_error=float(max(abs(s11-exact['s11']))),
                             measured_k=theta.real.tolist(), model_k=predicted_bloch.k.tolist(),
                             max_k_error=float(max(abs(theta.real-predicted_bloch.k))))
    save(); print('passive', json.dumps(result['passive']), flush=True)
    baseline = run_chain(runner, chain, root/'pump', [Tone(fp, power, ramp_time=args.ramp)], transient)
    actual = run_chain(runner, chain, root/'probe', [Tone(fp, power, ramp_time=args.ramp),
                        Tone(fs, probe, ramp_time=args.ramp)], transient)
    fitted = measure_chain(actual, chain, [fs, fi, fp, 3*fp], baseline)
    gain = float(abs(fitted['transmitted'][0]/fitted['incident'][0])**2)
    conversion = float(abs(fitted['transmitted'][1]/fitted['incident'][0])**2*fs/fi)
    model = four_wave_response(chain, fp, power, fs, probe, depleted=args.depleted,
                              harmonics=args.harmonics, sidebands=args.sidebands, quadrature=args.quadrature)
    pump_fit = measure_chain(actual, chain, [fp, 3*fp])
    model_pump = model.pump_transmission
    spectrum = measure_chain(actual, chain, model.frequency)
    result['active'] = dict(measured_gain=gain, model_gain=model.gain,
                            measured_conversion=conversion, model_conversion=model.idler_conversion,
                            gain_error_db=float(10*np.log10(gain/model.gain)),
                            conversion_error_db=float(10*np.log10(conversion/model.idler_conversion)),
                            model_max_phase=model.maximum_junction_phase, residual=model.residual,
                            harmonics=args.harmonics, sidebands=args.sidebands, quadrature=args.quadrature,
                            depleted=model.depleted, model_pump_transmission=model_pump,
                            measured_pump_transmission=float(abs(pump_fit['transmitted'][0]/pump_fit['incident'][0])**2),
                            third_harmonic_power_ratio=float(abs(pump_fit['transmitted'][1]/pump_fit['transmitted'][0])**2))
    result['active']['retained_frequency'] = model.frequency.tolist()
    result['active']['measured_output_power'] = (abs(spectrum['transmitted'])**2/(2*chain.port_impedance)).tolist()
    result['active']['model_output_power'] = model.transmitted_power.tolist()
    result['active']['pump_error_db'] = float(10*np.log10(result['active']['measured_pump_transmission']/model_pump))
    result['passed'] = bool(result['passive']['max_s21_error'] < .03 and
                            result['passive']['max_s11_error'] < .03 and
                            result['passive']['max_k_error'] < 1e-3 and
                            abs(result['active']['gain_error_db']) < .35 and
                            abs(result['active']['conversion_error_db']) < .5 and
                            abs(result['active']['pump_error_db']) < .1)
    save(); print('active', json.dumps(result['active']), flush=True)
    if not result['passed']:
        raise SystemExit('4WM circuit/model acceptance failed; inspect results before claiming validation.')


if __name__ == '__main__':
    main()
