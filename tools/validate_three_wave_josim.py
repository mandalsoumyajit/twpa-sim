"""Physical rf-SQUID/SNAIL 3WM circuit validation with JoSIM."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
from twpa.backends import JoSIM
from twpa.three_wave import (example_chain, chain_netlist, linear_response, bloch,
                             three_wave_response, quantum_noise)
from twpa.jpa_circuit import Tone, Transient
from twpa.josim_jtwpa import measure_chain, measure_bloch


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kind', choices=['rf-squid', 'snail'], required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--executable', default='josim-cli')
    p.add_argument('--wsl')
    p.add_argument('--corrector', type=int, default=0)
    p.add_argument('--cells', type=int, default=64)
    p.add_argument('--step', type=float, default=2e-12)
    p.add_argument('--stop', type=float, default=400e-9)
    p.add_argument('--start', type=float, default=300e-9)
    p.add_argument('--ramp', type=float, default=30e-9)
    p.add_argument('--bias-ramp', type=float, default=100e-9)
    p.add_argument('--bias-shape', choices=['smooth', 'linear'], default='smooth')
    p.add_argument('--pump-ratio', type=float, default=.25)
    p.add_argument('--probe-ratio', type=float, default=1e-6)
    p.add_argument('--junction-resistance', type=float, default=1e9)
    p.add_argument('--ground-resistance', type=float)
    p.add_argument('--loop-inductance', type=float, default=1e-13)
    p.add_argument('--harmonics', type=int, default=3)
    p.add_argument('--sidebands', type=int, default=1)
    p.add_argument('--quadrature', type=int, default=24)
    p.add_argument('--temperature', type=float, default=0.)
    p.add_argument('--depleted', action='store_true')
    p.add_argument('--engineered', action='store_true')
    args = p.parse_args()
    root = args.output.resolve(); root.mkdir(parents=True, exist_ok=False)
    chain = example_chain(args.kind, args.cells, junction_resistance=args.junction_resistance,
                          ground_resistance=args.ground_resistance, dispersion_engineered=args.engineered)
    c = chain.linear_chain
    transient = Transient(args.step, args.stop, args.start, 10e-12)
    runner = JoSIM(args.executable, args.wsl, timeout=600, corrector=args.corrector)
    fp, fs, fi = 1e9, .45e9, .55e9
    power = .5*c.port_impedance*(args.pump_ratio*c.critical_current)**2
    probe = power*args.probe_ratio
    result = dict(kind=args.kind, chain=asdict(chain), transient=asdict(transient),
                  simulator=runner.probe(), pump_frequency=fp, signal_frequency=fs,
                  pump_power=power, signal_power=probe, controls=dict(harmonics=args.harmonics,
                  sidebands=args.sidebands, quadrature=args.quadrature),
                  netlist_controls=dict(bias_ramp=args.bias_ramp, bias_shape=args.bias_shape,
                     snail_loop_inductance=args.loop_inductance),
                  source_ramp=args.ramp, engineered=args.engineered)
    def save():
        (root/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    def run(name, tones):
        deck = chain_netlist(chain, tones, transient, bias_ramp=args.bias_ramp,
                             snail_loop_inductance=args.loop_inductance, bias_shape=args.bias_shape)
        return runner.run(deck, root/name, metadata=dict(chain=asdict(chain),
                          kind=args.kind, tones=[asdict(t) for t in tones], transient=asdict(transient)))
    frequencies = [fs, fp, fi]
    passive = run('passive', [Tone(f, power*1e-7, ramp_time=args.ramp) for f in frequencies])
    m = measure_chain(passive, chain, frequencies)
    exact = linear_response(chain, frequencies)
    s21, s11 = m['transmitted']/m['incident'], m['reflected']/m['incident']
    k = measure_bloch(passive, chain, frequencies).real
    result['passive'] = dict(s21_real=s21.real.tolist(), s21_imag=s21.imag.tolist(),
        max_s21_error=float(max(abs(s21-exact['s21']))), max_s11_error=float(max(abs(s11-exact['s11']))),
        max_k_error=float(max(abs(k-bloch(chain, frequencies).k))), measured_k=k.tolist())
    save(); print('passive', json.dumps(result['passive']), flush=True)
    base = run('pump', [Tone(fp, power, ramp_time=args.ramp)])
    actual = run('probe', [Tone(fp, power, ramp_time=args.ramp), Tone(fs, probe, ramp_time=args.ramp)])
    m = measure_chain(actual, chain, [fs, fi, fp, 2*fp, 3*fp], base)
    gain = float(abs(m['transmitted'][0]/m['incident'][0])**2)
    conversion = float(abs(m['transmitted'][1]/m['incident'][0])**2*fs/fi)
    r = three_wave_response(chain, fp, power, fs, probe, depleted=args.depleted, **result['controls'])
    pm = measure_chain(actual, chain, [fp, 2*fp, 3*fp])
    pump = float(abs(pm['transmitted'][0]/pm['incident'][0])**2)
    middle = max(1, chain.cells//2)
    label = f'P(BJ{middle})' if args.kind=='rf-squid' else f'P(BJS{middle})'
    dc = float(np.mean(actual.trace(label)))
    predicted_dc = float(chain.element.equilibrium+r.dc_phase[middle])
    result['active'] = dict(measured_gain=gain, model_gain=r.gain, measured_conversion=conversion,
        model_conversion=r.idler_conversion, gain_error_db=float(10*np.log10(gain/r.gain)),
        conversion_error_db=float(10*np.log10(conversion/r.idler_conversion)),
        measured_pump_transmission=pump, model_pump_transmission=r.pump_transmission,
        pump_error_db=float(10*np.log10(pump/r.pump_transmission)), depleted=args.depleted,
        measured_mean_phase=dc, model_mean_phase=predicted_dc, dc_phase_error=dc-predicted_dc,
        maximum_phase=r.maximum_phase, residual=r.residual,
        second_harmonic_power_ratio=float(abs(pm['transmitted'][1]/pm['transmitted'][0])**2))
    q = quantum_noise(chain, fp, power, fs, junction_temperature=args.temperature,
                      ground_temperature=args.temperature, **result['controls'])
    result['quantum'] = dict(gain=q.gain, added_noise=q.added_noise, quantum_limit=q.quantum_limit,
        excess_noise=q.excess_noise, commutator_error=q.commutator_error,
        anomalous_commutator_error=q.anomalous_commutator_error,
        minimum_uncertainty_eigenvalue=q.minimum_uncertainty_eigenvalue,
        internal_temperature=args.temperature)
    a = result['active']; b = result['passive']
    result['passed'] = bool(b['max_s21_error']<.03 and b['max_s11_error']<.03 and b['max_k_error']<1e-3
        and abs(a['gain_error_db'])<.15 and abs(a['conversion_error_db'])<.2
        and abs(a['pump_error_db'])<.1 and abs(a['dc_phase_error'])<.003 and q.excess_noise>-1e-8)
    save(); print('active', json.dumps(a), flush=True); print('quantum', json.dumps(result['quantum']), flush=True)
    if not result['passed']:
        raise SystemExit('3WM circuit/model comparison failed; inspect data before claiming validation')


if __name__ == '__main__':
    main()
