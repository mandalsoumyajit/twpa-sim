"""Run the IDC electrostatics -> lossy circuit -> amplifier-noise pipeline."""
from twpa import preset, IDCGeometry, oxide_design, TLS, coupled

design, bias = preset('v2-tls')
coated, bare, fields = oxide_design(design, IDCGeometry(mesh=5e-9), oxide_thickness=30e-9)
print('Oxide participation:', fields.participation)
print('Capacitance scaling:', fields.energy_per_length/bare.energy_per_length)
# Illustrative material parameters, not measured device characterization.
tls = TLS(loss_tangent=.001, temperature=.02, critical_field=30., field_per_amp=2e8)
result = coupled(coated, bias, [3.5e9, 4e9, 4.5e9], tls)
print('Pump GHz:', result.pump_frequency/1e9)
print('Gain dB:', result.gain_db)
print('Added noise, quanta:', result.added_noise)
