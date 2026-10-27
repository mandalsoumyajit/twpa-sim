"""Build a KIT ladder, calculate dispersion, and compare weak/finite probes.

Run from the repository after `python -m pip install -e .`.
These electrical values are illustrative; no fabrication process is assumed.
"""
import json
import numpy as np
from twpa import Cell, Section, Design, Bias, linear, small_signal, nonlinear


def main():
    # Cell.capacitance is PER FINGER. Two 10 fF fingers give 20 fF per cell.
    # sqrt(50 pH / 20 fF) = 50 ohm in the low-frequency limit.
    cell = Cell(inductance=50e-12, capacitance=10e-15)
    design = Design((Section(cell, 1),), repeats=12000, name='howto-kit')

    # Currents are peak amperes, not RMS. DC bias enables the KIT 3WM term.
    bias = Bias(scaling_current=1e-3, dc_current=.2e-3, pump_current=20e-6)
    grid = np.linspace(0., 10e9, 10001)  # Include pump, signal and idler.
    passive = linear(design, grid, bias.impedance)

    # This uniform line has no designed stopband: select the pump explicitly.
    # For 3WM the idler is fp-fs, so the probe below has fi=4.2 GHz.
    pump, signal = 7e9, 2.8e9
    weak = small_signal(design, bias, [signal], frequency=grid, pump_frequency=pump)
    finite = nonlinear(design, bias, signal, bias.pump_current/10,
                       frequency=grid, pump_frequency=pump, samples=101)
    print(json.dumps(dict(family='KIT', cells=design.total_cells,
        signal_GHz=signal/1e9, idler_GHz=(pump-signal)/1e9,
        weak_gain_dB=float(weak.gain_db[0]), finite_gain_dB=float(10*np.log10(finite.gain)),
        output_pump_peak_A=float(abs(finite.currents[0,-1])),
        passive_max_power_balance_error=float(np.max(abs(abs(passive.s11)**2+abs(passive.s21)**2-1)))), indent=2))


if __name__ == '__main__':
    main()
