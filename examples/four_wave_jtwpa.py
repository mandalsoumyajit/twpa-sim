"""Compare propagation and amplification with and without resonant loading.

A finite line contains 384 Josephson junctions with node-to-ground capacitors.
The second circuit adds capacitively coupled inductor-capacitor resonators,
which alter dispersion and hence how contributions from successive cells
add in phase. Both use a 1 GHz pump and a 0.9 GHz signal; the principal
four-wave-mixing idler is 1.1 GHz.

PHIBAR is the reduced superconducting flux quantum. Setting critical current
to PHIBAR/(2.5 nH) chooses a zero-phase junction inductance of 2.5 nH.
The port resistance is matched to each circuit's Bloch impedance at the pump;
P=R*I_peak^2/2 sets available input power, not every branch's actual current.
See the printed labels for gain/conversion and compare complete solutions,
not just component counts. Start with examples/howto/05_four_wave_jtwpa.py.
Physical background: O'Brien et al.,
https://doi.org/10.1103/PhysRevLett.113.157001; docs/four-wave-jtwpa.md.
"""
from dataclasses import replace
import numpy as np
from twpa import jtwpa
from twpa.jpa_circuit import PHIBAR

base = jtwpa.Chain(critical_current=PHIBAR/2.5e-9, ground_capacitance=1e-12,
                   junction_capacitance=1e-12, cells=384)
load = jtwpa.ResonantLoad(coupling_capacitance=10e-15, capacitance=1e-12,
                          inductance=1/((2*np.pi*1.008e9)**2*1e-12))
for name, circuit in [('uniform', base), ('phase-matched', replace(base, resonator=load))]:
    circuit = replace(circuit, port_impedance=float(jtwpa.bloch(circuit, [1e9]).impedance[0].real))
    power = .5*circuit.port_impedance*(.4*circuit.critical_current)**2
    response = jtwpa.four_wave_response(circuit, 1e9, power, .9e9, power*1e-6)
    print(name, 'gain dB:', 10*np.log10(response.gain),
          'idler conversion:', response.idler_conversion,
          'junction phase bound:', response.maximum_junction_phase)
