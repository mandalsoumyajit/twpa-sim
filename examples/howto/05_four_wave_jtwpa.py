"""A finite 4WM junction chain, coherent response and quantum internal loss."""
import json
import numpy as np
from twpa import jtwpa


def main():
    # Ic is the physical single-junction critical current. Cg is per interior
    # node and Cj per junction. Both end shunts have half the interior value.
    chain = jtwpa.Chain(critical_current=1.3e-7, ground_capacitance=1e-12,
        junction_capacitance=1e-12, cells=64, port_impedance=50.,
        junction_resistance=100000., ground_resistance=100000.)
    fp, fs, pump_power = 1e9, .9e9, 7e-14
    # Four-wave mixing creates fi=2*fp-fs=1.1 GHz.
    response = jtwpa.four_wave_response(chain, fp, pump_power, fs,
                                        signal_power=pump_power*1e-6)
    # Temperatures in kelvin set bath occupations, without changing R or Ic.
    noise = jtwpa.quantum_noise(chain, fp, pump_power, fs,
                               junction_temperature=.02, ground_temperature=.02)
    print(json.dumps(dict(family='4WM JTWPA', gain_dB=float(10*np.log10(response.gain)),
        idler_photon_conversion=response.idler_conversion,
        added_quanta=noise.added_noise, finite_gain_limit=noise.quantum_limit,
        commutator_error=noise.commutator_error,
        maximum_phase_rad=response.maximum_junction_phase), indent=2))
    # Add jtwpa.ResonantLoad for resonant phase matching. For a stronger probe,
    # use depleted=True and increase sidebands/harmonics until gain converges.


if __name__ == '__main__':
    main()
