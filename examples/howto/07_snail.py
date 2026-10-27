"""A SNAIL chain with its full asymmetric potential and small-signal noise."""
import json
from twpa import three_wave as tw


def main():
    # The example has n=3 large JJs, a small/large Ic ratio alpha=0.3,
    # and flux=0.35 Phi0. Its reduced potential neglects geometric loop L.
    chain = tw.example_chain('snail', cells=64, dispersion_engineered=True)
    pump = .5*chain.port_impedance*(.18*chain.linear_chain.critical_current)**2
    result = tw.three_wave_response(chain, 1e9, pump, .45e9, pump*1e-6)
    # 3WM has fi=fp-fs=0.55 GHz; the full potential also retains residual Kerr.
    noise = tw.quantum_noise(chain, 1e9, pump, .45e9, junction_temperature=.02)
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
