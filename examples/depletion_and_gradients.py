"""Explore how assumed spatial inductance variation changes finite-signal gain.

The coupled-mode equations evolve coherent pump, signal and idler currents
along a kinetic-inductance line. A finite signal can deplete the pump. Here a
2% fractional inductance variation is applied through several normalized
spatial profiles, and the resulting forward power gain and remaining peak
pump current (A) are compared.

Names such as 'folded' and 'spiral' label mathematical variation profiles.
They do not insert bends, layout parasitics or a package field simulation.
This is a classical propagation study, not a quantum-noise calculation.
Start with examples/howto/01_kit.py; see docs/kit-twpa.md and its Malnou/Eom
references for physical context, and docs/numerical-methods.md for refinement.
"""
from twpa import preset, nonlinear, route_profile

design, bias = preset('v2-tls')
for mode, cycles in [('ideal', 1), ('linear', 1), ('sine', 7), ('folded', 1), ('folded', 2), ('spiral', 7)]:
    r = nonlinear(design, bias, 4.6e9, bias.pump_current/100,
                  gradient=.02, profile=route_profile(mode, cycles))
    print(mode, cycles, 'gain:', r.gain, 'output pump A:', abs(r.currents[0, -1]))
