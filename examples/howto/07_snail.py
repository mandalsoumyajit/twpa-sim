"""Calculate three-wave gain and noise in a line of asymmetric junction loops.

A SNAIL (superconducting nonlinear asymmetric inductive element) has one
small junction on one loop branch and several larger series junctions on
the other. Magnetic flux changes the equilibrium and nonlinear potential.
This example uses three large junctions, a small/large critical-current ratio
of 0.3 and applied flux equal to 0.35 flux quanta in a synthetic 64-cell line.
The 1 GHz pump couples a 0.45 GHz signal to a 0.55 GHz idler.

The output gives equilibrium small-junction phase (radians), normalized
potential derivatives c2/c3/c4 (linear curvature, cubic and quartic terms),
transmitted power gain and idler photon conversion. Added noise is compared
with the gain-dependent quantum bound. A small commutator defect checks
quantum algebra, not physical stability or convergence of omitted harmonics.

Read HOWTO.md and docs/three-wave-josephson.md. Element background:
Frattini et al., https://doi.org/10.1063/1.4984142; traveling-wave context:
Roudsari et al., https://doi.org/10.1063/5.0127690.
"""
import json
from twpa import three_wave as tw


def main():
    # The example has n=3 large JJs, a small/large Ic ratio alpha=0.3,
    # and flux=0.35 Phi0. Its reduced potential neglects geometric loop L.
    chain = tw.example_chain('snail', cells=64, dispersion_engineered=True)
    # I_c*c2 sets the small-signal curvature scale used to choose this drive;
    # it is not the physical large junction's critical current. Convert an
    # incoming peak-current scale to available power with P=R*I_peak^2/2.
    pump_power = .5*chain.port_impedance*(.18*chain.linear_chain.critical_current)**2
    result = tw.three_wave_response(chain, 1e9, pump_power, .45e9, pump_power*1e-6)
    # 3WM has fi=fp-fs=0.55 GHz; the full potential also retains residual Kerr.
    noise = tw.quantum_noise(chain, 1e9, pump_power, .45e9, junction_temperature=.02)
    print(json.dumps(dict(family='SNAIL 3WM',
        equilibrium_small_JJ_phase_rad=chain.element.equilibrium,
        potential_c2_c3_c4=tw.coefficients(chain.element).tolist(), gain=result.gain,
        idler_photon_conversion=result.idler_conversion,
        added_quanta=noise.added_noise, finite_gain_limit=noise.quantum_limit,
        commutator_error=noise.commutator_error), indent=2))
    # Start a custom device with SNAIL(Ic_large, alpha, n, Cj_large, ...).
    # All JJs in a cell share its junction bath temperature. Quantum noise is
    # about the pump-only state; classical depleted=True is a separate option.


if __name__ == '__main__':
    main()
