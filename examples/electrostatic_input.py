"""Turn a small electrostatic cross-section into an illustrative cell value.

The grid holds a dielectric between boundaries fixed at 0 V and 1 V.
Solving its potential gives electric energy and capacitance per invariant
length. A selected strip of cells marks a region whose energy participation
is measured. This region label does not independently change permittivity.
Multiplying capacitance per length by an assumed finger length then creates
a toy electrical capacitance; it is not a full interdigitated-capacitor or
amplifier-layout extraction.

er means dimensionless relative permittivity. fixed contains prescribed
voltages (NaN means unknown), and oxide is a Boolean participation mask.
dx/dy are mesh spacings in metres. The printed Cell has inductance in H,
capacitance in F and dimensionless participation. See docs/kit-twpa.md and
examples/howto/02_tls_loss.py for the dielectric-loss pipeline and references.
"""
import numpy as np
from twpa import electrostatic, Cell

er = np.ones((41, 16))*4
fixed = np.full(er.shape, np.nan)
fixed[0, :] = 0
fixed[-1, :] = 1
oxide = np.zeros_like(er, bool)
oxide[1:4, :] = True
result = electrostatic(er, fixed, oxide, dx=2.5e-9, dy=2.5e-9)
# For this illustrative invariant cross-section only:
finger_length = 72e-6
cell = Cell(54e-12, result.capacitance_per_length*finger_length,
            stub_inductance=1.3e-9, participation=result.participation)
print(cell)
