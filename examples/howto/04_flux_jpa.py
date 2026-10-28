"""Amplify one phase component and squeeze the other in a flux-pumped JPA.

An effective magnetic-flux pump at 10 GHz modulates a 5 GHz resonator.
Signal and idler coincide at 5 GHz (degenerate operation), so the useful
observables are the cosine-like and sine-like field quadratures relative
to the pump. One combination is amplified and the orthogonal one reduced.

The output lists their two power gains and symmetrized spectral noise
variances. Compare noise with vacuum=0.5; below 0.5 means quadrature
squeezing in this normalization. A single phase-preserving added-noise
figure is undefined here. The effective pair rate is assumed, not extracted
from a physical coil's power or flux.

See docs/resonant-jpa.md and its Eichler/Wallraff and Caves references:
https://doi.org/10.1140/epjqt2 and https://doi.org/10.1103/PhysRevD.26.1817.
"""
import json
import numpy as np
from twpa import jpa


def main():
    r = jpa.Resonator(frequency=5e9, kappa_external=2*np.pi*10e6,
                      kappa_internal=2*np.pi*.25e6, internal_temperature=.02)
    # Physical pump frequency is twice the frame/signal frequency.
    # pair_drive is an effective angular rate, NOT a pump power in watts.
    response = jpa.flux_response(r, pump_frequency=10e9,
        pair_drive=.3*r.kappa, signal_frequency=[5e9])
    # At degeneracy there is one mode, and the two quadratures have different
    # gains. response.added_noise is intentionally NaN in this case.
    q = jpa.quadratures(response)
    print(json.dumps(dict(family='Effective flux JPA',
        quadrature_power_gains=q.power_gains.tolist(),
        output_noise_variances=q.noise_variances.tolist(),
        vacuum_variance=.5), indent=2))
    # To map a physical SQUID flux drive to this effective rate, use the
    # SQUIDJPA adapter/JoSIM comparison described in docs/josim-backend.md.


if __name__ == '__main__':
    main()
