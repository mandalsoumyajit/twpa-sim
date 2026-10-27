"""Simple cross-section example, not a KITWPA geometry extraction."""
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
