"""Connect IDC electrostatics to the coupled KIT/TLS gain and noise solver."""
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

    # These TLS numbers are illustrative. Calibrate field_per_amp and the
    # material loss tangent for your device; linewidth is phenomenological.
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
