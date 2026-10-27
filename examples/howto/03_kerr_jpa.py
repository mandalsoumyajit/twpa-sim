"""Select a stable Kerr-JPA pump branch and calculate gain and loss noise."""
import json
import numpy as np
from twpa import jpa


def main():
    # Frequencies are Hz. Linewidths and Kerr coefficient are ANGULAR rates.
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
