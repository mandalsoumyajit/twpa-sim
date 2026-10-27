"""Compare historical electrical presets. Run after pip install -e ."""
import json
import numpy as np
from twpa import preset, small_signal

for name in ('v1-target', 'v1-cpw', 'v1-idc', 'v2'):
    design, bias = preset(name)
    result = small_signal(design, bias, np.linspace(1.5e9, 8e9, 400))
    print(json.dumps(dict(design=name, cells=design.total_cells,
                          pump_ghz=result.pump_frequency/1e9,
                          peak_gain_db=float(np.nanmax(result.gain_db)))))
