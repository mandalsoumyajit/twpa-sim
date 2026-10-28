"""Use magnetic bias to enable three-wave mixing in a line of one-JJ loops.

An rf-SQUID (radio-frequency superconducting quantum interference device)
has one Josephson junction and a geometric loop inductance. Applied flux
shifts its equilibrium phase, giving the asymmetric response needed for
three-wave mixing: a 1 GHz pump generates signal/idler response at 0.45 and
0.55 GHz. The loop's own circulating current screens the applied flux.

We solve a synthetic 64-cell line, then calculate quantum noise including
junction and ground-resistor baths at 20 mK. Gain is a transmitted power
ratio. flux_Phi0 is applied flux divided by the flux quantum. The potential
coefficients c2/c3/c4 are dimensionless derivatives at equilibrium: linear
curvature, cubic nonlinearity, and quartic nonlinearity. The printed DC shift
is pump-induced static phase rectification in radians.

See HOWTO.md and docs/three-wave-josephson.md. Physical background:
Zorin, https://doi.org/10.1103/PhysRevApplied.6.034006.
"""
import json
from twpa import three_wave as tw


def main():
    # beta_L = L_loop*I_c/phi_bar = 0.8, where L_loop is loop inductance,
    # I_c is the junction critical current, and phi_bar is the flux quantum
    # divided by 2*pi. Below one, the model supports a unique static branch.
    # The engineered loading suppresses unwanted sum-frequency propagation.
    chain = tw.example_chain('rf-squid', cells=64, dispersion_engineered=True,
                              junction_resistance=20000., ground_resistance=100000.)
    c2, c3, c4 = tw.coefficients(chain.element)
    # The normalized potential for phase displacement x is
    # U/E_J = constant+c2*x^2/2+c3*x^3/6+c4*x^4/24+... .
    # U is energy and E_J=phi_bar*I_c; c3 permits three-wave conversion.
    # linear_chain.critical_current is Ic*c2, a CURVATURE scale. It is not
    # the physical rf-SQUID junction's switching current.
    available_pump = .5*chain.port_impedance*(.18*chain.linear_chain.critical_current)**2
    # P=R*I_peak^2/2 converts a reference incoming-wave current to available
    # power in watts. Reflections mean it does not prescribe every JJ current.
    result = tw.three_wave_response(chain, 1e9, available_pump, .45e9,
                                    signal_power=available_pump*1e-6)
    noise = tw.quantum_noise(chain, 1e9, available_pump, .45e9,
                             junction_temperature=.02, ground_temperature=.02)
    print(json.dumps(dict(family='rf-SQUID 3WM', flux_Phi0=chain.element.flux_bias,
        potential_c2_c3_c4=[c2,c3,c4], gain=result.gain,
        idler_photon_conversion=result.idler_conversion,
        maximum_DC_shift_rad=float(max(abs(result.dc_phase))),
        added_quanta=noise.added_noise, finite_gain_limit=noise.quantum_limit), indent=2))
    # Configure tw.RFSQUID(Ic, Lloop, Cj, ..., flux_bias=...) and ThreeWaveChain
    # directly to replace the example with your own physical parameters.


if __name__ == '__main__':
    main()
