"""Find the reflected gain of a pumped nonlinear microwave resonator.

A Josephson parametric amplifier (JPA) is represented here by one effective
mode, rather than individual junctions. Its Kerr coefficient makes resonance
shift with stored photon number. We first find steady pump states, choose
a stable one, and calculate the response to a weak signal plus its idler.

The cavity photon count describes pump buildup. Signal gain is reflected
power gain in dB; added quanta are symmetrized noise referred to the input.
The 0.5 MHz signal offset gives an idler 0.5 MHz below the pump. All values
are illustrative; internal temperature sets the loss bath's fluctuations.

Read docs/resonant-jpa.md and HOWTO.md. Physical treatment:
Eichler and Wallraff, EPJ Quantum Technology 1, 2 (2014),
https://doi.org/10.1140/epjqt2.
"""
import json
import numpy as np
from twpa import jpa


def main():
    # The mode is at 5 GHz. kappa_external and kappa_internal are energy
    # decay rates into the port and internal bath. Convert a linewidth in Hz
    # to an angular rate with 2*pi. Field amplitude decays at total kappa/2.
    # Negative Kerr shifts angular resonance downward by 2*pi*10 kHz per
    # additional cavity photon. internal_temperature is 0.020 kelvin.
    resonator = jpa.Resonator(frequency=5e9, kappa_external=2*np.pi*10e6,
        kappa_internal=2*np.pi*.25e6, kerr=-2*np.pi*10e3,
        internal_temperature=.02)
    roots = jpa.pump_points(resonator, frequency=4.994e9, power=2e-13)
    # 'lower' chooses the lowest-occupation STABLE root. Do not select a
    # bistable branch merely by taking the first polynomial root.
    pump = jpa.select_pump_point(roots, branch='lower')
    signal = pump.frequency+.5e6
    response = jpa.kerr_response(resonator, pump, [signal])
    # In Kerr pumping fi=2*fp-fs. Stay away from fs=fp for this noise figure.
    print(json.dumps(dict(family='Kerr JPA', pump_roots=len(roots),
        selected_branch=pump.stability, cavity_photons=pump.photons,
        signal_GHz=signal/1e9, idler_GHz=float(response.idler_frequency[0]/1e9),
        gain_dB=float(response.gain_db[0]), added_quanta=float(response.added_noise[0])), indent=2))


if __name__ == '__main__':
    main()
