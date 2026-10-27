"""Coupled TLS/amplifier example; loss parameters are illustrative."""
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
