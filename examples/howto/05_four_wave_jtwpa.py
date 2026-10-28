"""Follow signal amplification and loss noise through a finite junction chain.

There are 64 Josephson junctions in series and capacitors from the nodes to
ground. Junction capacitance and loss resistance are retained too. A 1 GHz
pump couples a 0.9 GHz signal to a 1.1 GHz idler by four-wave mixing (4WM).
Physical source/load resistors allow reflected waves in the circuit solution.

The classical call finds transmitted weak-signal gain and generated idler.
The quantum call propagates small fluctuations from ports and internal loss
baths about the pump. Compare added_quanta with finite_gain_limit; the latter
depends on gain and approaches 0.5 only at high gain. Idler conversion counts
photons, not watts. commutator_error is a quantum consistency defect, while
maximum_phase_rad is a phase-excursion diagnostic, not a switching threshold.

This short illustrative line is not optimized for large gain. See HOWTO.md,
docs/four-wave-jtwpa.md and docs/jtwpa-quantum-noise.md. Background:
O'Brien et al., https://doi.org/10.1103/PhysRevLett.113.157001;
Caves, https://doi.org/10.1103/PhysRevD.26.1817.
"""
import json
import numpy as np
from twpa import jtwpa


def main():
    # Critical current I_c sets junction inductance. Ground capacitance C_g
    # connects each interior node to ground, while junction capacitance C_j
    # lies across each junction. End ground shunts have half the interior
    # value. The two finite 100 kohm resistances are dissipative loss paths.
    chain = jtwpa.Chain(critical_current=1.3e-7, ground_capacitance=1e-12,
        junction_capacitance=1e-12, cells=64, port_impedance=50.,
        junction_resistance=100000., ground_resistance=100000.)
    pump_frequency, signal_frequency, pump_power = 1e9, .9e9, 7e-14
    # Frequencies are Hz, available incident power is W. The signal is a
    # million times weaker in power, so this is a small-signal calculation.
    # Four-wave mixing: idler frequency = 2*pump_frequency-signal_frequency.
    response = jtwpa.four_wave_response(chain, pump_frequency, pump_power, signal_frequency,
                                        signal_power=pump_power*1e-6)
    # Temperatures in kelvin set bath occupations, without changing R or Ic.
    noise = jtwpa.quantum_noise(chain, pump_frequency, pump_power, signal_frequency,
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
