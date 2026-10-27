"""Flux-biased rf-SQUID 3WM, including static screening and resistor baths."""
import json
from twpa import three_wave as tw


def main():
    # This synthetic example uses beta_L=0.8 (<1: unique supported DC branch).
    # The engineered loading suppresses unwanted sum-frequency propagation.
    chain = tw.example_chain('rf-squid', cells=64, dispersion_engineered=True,
                              junction_resistance=20000., ground_resistance=100000.)
    c2, c3, c4 = tw.coefficients(chain.element)
    # linear_chain.critical_current is Ic*c2, a CURVATURE scale. It is not
    # the physical rf-SQUID junction's switching current.
    available_pump = .5*chain.port_impedance*(.18*chain.linear_chain.critical_current)**2
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
