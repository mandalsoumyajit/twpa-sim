"""Does a biased superconducting transmission line amplify a traveling signal?

The line's kinetic inductance depends on current. A direct-current bias
enables three-wave mixing: the 7 GHz pump couples a 2.8 GHz signal to a
4.2 GHz idler. We compare an undepleted weak-signal calculation with one
that evolves all three finite classical currents, allowing pump depletion.

The printed gains are transmission power gains in dB (0 dB is unity).
The remaining pump current is a peak value in amperes. The passive balance
error checks reflected plus transmitted power against incident power.
This illustrative electrical ladder is not an extracted mask or a calibrated
material model. No quantum noise is computed by the finite-current solve.

Run from the repository after ``python -m pip install -e .``.
Read HOWTO.md and docs/kit-twpa.md for the equations and output conventions.
Physical background: Malnou et al., PRX Quantum 2, 010302 (2021),
https://doi.org/10.1103/PRXQuantum.2.010302.
"""
import json
import numpy as np
from twpa import Cell, Section, Design, Bias, linear, small_signal, nonlinear


def main():
    # Cell.capacitance is PER FINGER. Two 10 fF fingers give 20 fF per cell.
    # sqrt(L/C) is the low-frequency wave impedance: the inductance L=50 pH
    # and total ground capacitance C=20 fF give a 50 ohm line.
    cell = Cell(inductance=50e-12, capacitance=10e-15)
    design = Design((Section(cell, 1),), repeats=12000, name='howto-kit')

    # Currents are peak amperes, not root-mean-square values. The scaling
    # current sets the strength of current-dependent inductance; it is not
    # a prediction of the film's measured switching current.
    bias = Bias(scaling_current=1e-3, dc_current=.2e-3, pump_current=20e-6)
    grid = np.linspace(0., 10e9, 10001)  # Include pump, signal and idler.
    passive = linear(design, grid, bias.impedance)

    # This uniform line has no designed stopband: select the pump explicitly.
    # For 3WM the idler is fp-fs, so the probe below has fi=4.2 GHz.
    pump, signal = 7e9, 2.8e9
    weak = small_signal(design, bias, [signal], frequency=grid, pump_frequency=pump)
    finite = nonlinear(design, bias, signal, bias.pump_current/10,
                       frequency=grid, pump_frequency=pump, samples=101)
    # Compare the weak and finite gains before interpreting a depletion
    # trend. Recheck propagation integration and frequency-grid resolution
    # when changing line length or driving closer to a band edge.
    print(json.dumps(dict(family='KIT', cells=design.total_cells,
        signal_GHz=signal/1e9, idler_GHz=(pump-signal)/1e9,
        weak_gain_dB=float(weak.gain_db[0]), finite_gain_dB=float(10*np.log10(finite.gain)),
        output_pump_peak_A=float(abs(finite.currents[0,-1])),
        passive_max_power_balance_error=float(np.max(abs(abs(passive.s11)**2+abs(passive.s21)**2-1)))), indent=2))


if __name__ == '__main__':
    main()
