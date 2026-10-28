"""Connect surface-oxide electric fields to KIT gain and internal-loss noise.

An interdigitated capacitor (IDC) stores some of its electric energy in an
oxide that can contain absorbing two-level systems (TLS). A coated/uncoated
edge-domain calculation supplies participation and a capacitance ratio. The
coupled propagation model uses those values to evolve the attenuating pump,
weak signal/idler conversion and thermal/vacuum bath noise together.

Participation and capacitance scaling are dimensionless. Pump frequency is
GHz, gain is in dB and added noise is in symmetrized quanta referred to the
signal input. The field/current mapping and material loss parameters need
calibration; this small domain is not an extraction of an entire chip.
Begin with examples/howto/02_tls_loss.py. See docs/kit-twpa.md and its
Martinis/Gao references for dielectric-loss context.
"""
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
