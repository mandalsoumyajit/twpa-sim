"""How does a lossy surface oxide affect a kinetic-inductance amplifier?

An interdigitated capacitor (IDC) has interleaved metal fingers. Its edge
electric field stores some energy in oxide, where two-level systems (TLS)
can absorb energy. We solve a small coated/uncoated edge domain, use its
energy ratio to adjust capacitance, then propagate pump loss, weak-signal
conversion and bath noise together in the amplifier.

Participation is the oxide's fraction of modeled electric energy, not a
full-chip measurement. Capacitance scale is a coated/bare energy ratio.
Gain is in dB; added quanta are noise referred to the signal input after
subtracting that input's original noise. Material and field/current conversion
values here are assumptions to calibrate, not measured process parameters.

See docs/kit-twpa.md and HOWTO.md. Background: Martinis et al.,
https://doi.org/10.1103/PhysRevLett.95.210503, and Gao et al.,
https://doi.org/10.1063/1.2906373. These do not calibrate this example.
"""
from dataclasses import replace
import json
import numpy as np
from twpa import preset, IDCGeometry, oxide_design, TLS, coupled


def main():
    # A historical electrical preset is a convenient starting point, not a
    # universal material model. Shorten it here to keep the tutorial quick.
    design, bias = preset('v2-tls')
    design = replace(design, repeats=200)
    # The small half-gap electrostatic domain needs mesh AND domain studies
    # before using it quantitatively. The 5 nm mesh is an introductory run.
    geometry = IDCGeometry(mesh=5e-9)
    coated, bare, fields = oxide_design(design, geometry, oxide_thickness=30e-9)

    # loss_tangent is dimensionless dielectric dissipation, temperature is K,
    # critical_field is the saturation electric field in V/m, and
    # field_per_amp maps line current (A) to local field (V/m).
    # linewidth is the phenomenological pump-saturation hole width in Hz.
    tls = TLS(loss_tangent=.001, temperature=.02, critical_field=30.,
              field_per_amp=2e8, linewidth=1e6)
    grid = np.linspace(0., 11.5e9, 16385)
    # The solver selects a supported pump near a stopband in this preset.
    # Pump attenuation, local signal/idler loss and bath noise evolve together.
    result = coupled(coated, bias, [4e9], tls, frequency=grid)
    print(json.dumps(dict(family='KIT + TLS', participation=fields.participation,
        capacitance_scale=fields.energy_per_length/bare.energy_per_length,
        pump_GHz=result.pump_frequency/1e9, gain_dB=float(result.gain_db[0]),
        added_quanta=float(result.added_noise[0]),
        output_pump_peak_A=float(result.pump_output_current[0])), indent=2))


if __name__ == '__main__':
    main()
