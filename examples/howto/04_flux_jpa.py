"""A flux-pumped resonator at degeneracy: use quadratures, not two vacua."""
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
