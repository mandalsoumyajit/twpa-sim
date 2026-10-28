"""Calculate KIT gain and noise with a prescribed lossy-dielectric fraction.

This version starts with electrical cells and assigns 3% electric-energy
participation to a lossy dielectric rather than solving its geometry. A
two-level-system (TLS) loss model supplies field-dependent absorption. The
coupled solver evolves pump attenuation, weak signal/idler conversion and
the loss environment's fluctuations together.

The output gives pump/signal frequencies in GHz, gain in dB, and added noise
in symmetrized quanta referred to the input. Material loss tangent,
temperature, saturation field and current-to-field conversion are illustrative
assumptions. They are not determined by assigning participation=0.03.
Start with examples/howto/02_tls_loss.py; docs/kit-twpa.md explains the
pipeline and connects it to the Martinis/Gao dielectric-loss references.
"""
from dataclasses import replace
import json
from twpa import preset, TLS, coupled

design, bias = preset('v2')
design = replace(design, sections=tuple(
    replace(s, cell=replace(s.cell, participation=.03)) for s in design.sections))
tls = TLS(loss_tangent=.001, temperature=.02, critical_field=10., field_per_amp=4e8)
result = coupled(design, bias, [3e9, 4e9, 5e9], tls)
print(json.dumps(dict(pump_ghz=result.pump_frequency/1e9,
                      signal_ghz=(result.signal_frequency/1e9).tolist(),
                      gain_db=result.gain_db.tolist(),
                      added_noise_quanta=result.added_noise.tolist(),
                      pump_output_amp=result.pump_output_current.tolist()), indent=2))
