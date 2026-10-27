"""Parameterized full nonlinear CMEs and spatial phase-integral gradients."""
from twpa import preset, nonlinear, route_profile

design, bias = preset('v2-tls')
for mode, cycles in [('ideal', 1), ('linear', 1), ('sine', 7), ('folded', 1), ('folded', 2), ('spiral', 7)]:
    r = nonlinear(design, bias, 4.6e9, bias.pump_current/100,
                  gradient=.02, profile=route_profile(mode, cycles))
    print(mode, cycles, 'gain:', r.gain, 'output pump A:', abs(r.currents[0, -1]))
