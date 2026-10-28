"""Compare weak-signal gain across historical kinetic-inductance line models.

Each preset supplies electrical cells and a current bias. For each line we
choose the model's supported pump near a loading stopband and sweep signal
frequency, then print the largest forward power gain in decibels. The named
presets represent different electrical assumptions and lengths; they are not
automatically equivalent mask designs or a ranking of fabricated amplifiers.

The pump frequency is in GHz and cell count is the repeated electrical length.
Inspect the full dispersion/gain spectrum before treating a peak as useful
bandwidth. See examples/howto/01_kit.py, docs/kit-twpa.md and
docs/source-audit.md for physical references and parameter provenance.
Run from the repository after ``python -m pip install -e .``.
"""
import json
import numpy as np
from twpa import preset, small_signal

for name in ('v1-target', 'v1-cpw', 'v1-idc', 'v2'):
    design, bias = preset(name)
    result = small_signal(design, bias, np.linspace(1.5e9, 8e9, 400))
    print(json.dumps(dict(design=name, cells=design.total_cells,
                          pump_ghz=result.pump_frequency/1e9,
                          peak_gain_db=float(np.nanmax(result.gain_db)))))
