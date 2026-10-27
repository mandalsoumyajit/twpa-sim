"""Finite-chain 4WM examples; no external simulator required."""
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
