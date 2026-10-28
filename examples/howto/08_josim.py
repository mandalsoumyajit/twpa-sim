"""Generate a physical JoSIM deck; run it only when an executable is supplied.

The earlier tutorials solve amplifier equations in Python. Here we translate
a circuit into a netlist (a text list of components and connections) for
JoSIM, an independent classical time-domain simulator. Its junctions use
the resistively and capacitively shunted junction (RCSJ) model: sine-law
supercurrent, junction capacitance and resistive damping. Flux-biased loops
include physical inductors and coils.

The drive is a weak classical probe, not a pump/probe gain measurement.
An executed run produces voltage and phase traces plus its exact input and
provenance. Pumped validation requires the separate drivers' pump-only
subtraction, port-wave extraction and timestep/settling studies.
Classical traces do not independently measure quantum vacuum noise.

Default execution is a dry run and needs no simulator. Use a fresh output
directory for each actual simulation; see HOWTO.md for native/WSL commands.
WSL means Windows Subsystem for Linux. See docs/josim-backend.md and
Delport et al., https://doi.org/10.1109/TASC.2019.2897312.
"""
import argparse
from pathlib import Path
from twpa import jtwpa, three_wave
from twpa.backends import JoSIM
from twpa.jpa_circuit import PHIBAR, ParallelJPA, Tone, Transient, jpa_netlist
import numpy as np


def make_deck(family):
    # Integrate with a 1 ps timestep; save only after 300 ns of settling and
    # stop at 400 ns. output_step is the saved-trace spacing, not the internal
    # integration timestep. These must be refined for a quantitative study.
    transient = Transient(step=1e-12, stop=400e-9, start=300e-9, output_step=10e-12)
    tones = [Tone(1e9, 1e-22, ramp_time=30e-9)]  # Weak classical probe.
    if family == 'jpa':
        # Choose a toy 1 GHz parallel resonator with total capacitance 1 nF.
        # Its linear Josephson inductance is L_J=1/(C_total*omega_0^2),
        # where omega_0=2*pi*1 GHz. The required critical current is
        # I_c=phi_bar/L_J, and PHIBAR is the reduced flux quantum phi_bar.
        # The constructor supplies 1 pF junction C in addition to 999 pF shunt C.
        circuit = ParallelJPA(PHIBAR*1e-9*(2*np.pi*1e9)**2, 999e-12)
        return jpa_netlist(circuit, tones, transient)
    if family == 'jtwpa':
        circuit = jtwpa.Chain(1.3e-7, 1e-12, 1e-12, cells=4)
        return jtwpa.chain_netlist(circuit, tones, transient)
    # Real junction loops and mutual-inductance flux bias, not a fictitious
    # flux input attached directly to a Josephson-junction model.
    circuit = three_wave.example_chain(family, cells=4, dispersion_engineered=True)
    return three_wave.chain_netlist(circuit, tones, transient)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family', choices=['jpa','jtwpa','rf-squid','snail'], default='snail')
    p.add_argument('--write-netlist', type=Path)
    p.add_argument('--executable')
    p.add_argument('--wsl')
    p.add_argument('--corrector', type=int, default=0)
    p.add_argument('--run-directory', type=Path)
    args = p.parse_args()
    deck = make_deck(args.family)
    if args.write_netlist:
        args.write_netlist.parent.mkdir(parents=True, exist_ok=True)
        # Preserve an existing deck rather than silently replacing it.
        with args.write_netlist.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(deck)
        print('Saved netlist:', args.write_netlist)
    if not args.executable:
        print(f'Dry run: {args.family}, {len(deck.splitlines())} netlist lines; no JoSIM execution.')
        return
    if args.run_directory is None:
        p.error('--run-directory is required with --executable')
    runner = JoSIM(args.executable, args.wsl, corrector=args.corrector)
    run = runner.run(deck, args.run_directory, metadata=dict(tutorial_family=args.family))
    print('Classical waveform:', run.directory/'output.csv')
    # Transients validate classical response. Quantum vacuum noise must be
    # checked using input-output identities, not classical JoSIM noise traces.


if __name__ == '__main__':
    main()
