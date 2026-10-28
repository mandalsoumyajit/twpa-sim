"""Explore reflected gain and squeezing in effective resonant JPA models.

A nonlinear microwave resonator stores pump energy and mixes a weaker signal
with an idler. The Kerr case pumps through the signal port; the flux case
uses an assumed pair-coupling rate for resonance modulation at twice the
signal frequency. Their measured outputs are reflection gain and, when
signal/idler coincide, two phase-quadrature gains and noise variances.

Frequencies are Hz, power is W, and kappa/Kerr/pair-drive quantities are
angular rates in rad/s. Cavity occupation counts stored pump photons; wave
amplitude squared counts photons/s. Noise uses vacuum variance 0.5.
See examples/howto/03_kerr_jpa.py and 04_flux_jpa.py before running this sweep.
Physical background: Eichler and Wallraff, https://doi.org/10.1140/epjqt2;
the guide docs/resonant-jpa.md defines all equation symbols and noise outputs.

Run: python examples/resonant_jpa.py [--plot gain.png]
All example parameters are illustrative, not a characterized device.
"""
import argparse
import json
import numpy as np
from scipy.constants import h
from twpa import jpa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot', help='Optional figure filename; requires matplotlib')
    args = parser.parse_args()
    resonator = jpa.Resonator(frequency=5e9, kappa_external=2*np.pi*10e6,
                             kappa_internal=2*np.pi*.25e6, kerr=-2*np.pi*10e3)
    # Choose a target pump occupation, then compute its incident power.
    photons = 300.
    delta = 2*resonator.kerr*photons
    fp = resonator.frequency+delta/(2*np.pi)
    power = photons*((resonator.kappa/2)**2+(delta-resonator.kerr*photons)**2)*h*fp/resonator.kappa_external
    point = jpa.select_pump_point(jpa.pump_points(resonator, fp, power))
    frequency = fp+np.linspace(-10e6, 10e6, 501)
    response = jpa.kerr_response(resonator, point, frequency)
    squeezing = jpa.quadratures(jpa.flux_response(
        resonator, 2*resonator.frequency, .3*resonator.kappa, [resonator.frequency]))
    print(json.dumps(dict(pump_frequency_ghz=fp/1e9, pump_power_dbm=10*np.log10(power/1e-3),
                          pump_photons=point.photons, stability=point.stability,
                          max_sideband_gain_db=float(np.max(response.gain_db[~response.degenerate])),
                          flux_quadrature_gains_db=(10*np.log10(squeezing.power_gains)).tolist(),
                          flux_quadrature_noise_db_relative_vacuum=(10*np.log10(2*squeezing.noise_variances)).tolist()), indent=2))
    if args.plot:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), layout='constrained')
        axes[0].plot((frequency-fp)/1e6, response.gain_db)
        axes[0].set(xlabel='Signal offset from pump (MHz)', ylabel='Sideband gain (dB)',
                    title='Kerr-pumped resonant JPA')
        valid = ~response.degenerate
        axes[1].plot((frequency[valid]-fp)/1e6, response.added_noise[valid])
        axes[1].set(xlabel='Signal offset from pump (MHz)', ylabel='Added noise (quanta)',
                    title='Vacuum input idler, internal loss included')
        for axis in axes:
            axis.grid(alpha=.3)
        fig.savefig(args.plot, dpi=180)


if __name__ == '__main__':
    main()
