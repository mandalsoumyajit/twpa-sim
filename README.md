# TWPA Sim

Python models for traveling-wave and resonant parametric amplifiers.

Model dispersion, scattering, nonlinear gain, pump depletion and dielectric
loss in periodic transmission lines, or calculate operating points, gain and
noise for a single-mode resonant Josephson parametric amplifier (JPA).
Circuit parameters and operating conditions are explicit inputs.

**Implemented:** periodic ladder circuits and kinetic-inductance TWPA presets;
3WM propagation with Kerr phase shifts; coupled TLS loss and amplifier noise;
Kerr-pumped resonant JPAs; effective flux-pumped resonant amplifiers;
4WM Josephson junction chains with optional resonant phase matching,
multiharmonic mixing, pump depletion and quantum noise with internal loss.
An optional JoSIM backend generates
physical RCSJ circuits and has been used to check both JPA and JTWPA models.

## Installation

Python 3.10+, NumPy and SciPy are required. The core package does not require
MATLAB, KLayout or JoSIM.

```powershell
git clone https://github.com/mandalsoumyajit/twpa-sim.git
cd twpa-sim
python -m venv .venv
.\.venv\Scripts\python -m pip install -e .
.\.venv\Scripts\python -m unittest discover -s tests -v
```

On macOS/Linux, use `.venv/bin/python` instead of `.\.venv\Scripts\python`.
Install `.[plot]` for the optional example plots.

## Resonant JPA

```python
import numpy as np
from twpa import jpa

r = jpa.Resonator(frequency=5e9, kappa_external=2*np.pi*10e6,
                  kappa_internal=2*np.pi*.25e6, kerr=-2*np.pi*10e3)
points = jpa.pump_points(r, frequency=4.994e9, power=2e-13)
p = jpa.select_pump_point(points, branch='lower')
response = jpa.kerr_response(r, p, p.frequency+np.array([-.5e6, .5e6]))
print(response.gain_db, response.added_noise)
```

All pump roots are reported with their stability. Bistable operation requires
an explicit branch choice. The flux-pumped model also exposes quadrature gains
and output squeezing at degeneracy. Frequencies use Hz; decay rates and Kerr
coefficients use rad/s; incident power uses watts.

See the [JPA guide](docs/resonant-jpa.md) for normalization, noise conventions,
classical probe transients and limitations. Run
`python examples/resonant_jpa.py`, optionally adding `--plot gain.png`.

## Four-wave Josephson TWPA

```python
import numpy as np
from twpa import jtwpa

chain = jtwpa.Chain(critical_current=1.3e-7, ground_capacitance=1e-12,
                    junction_capacitance=1e-12, cells=256)
r = jtwpa.four_wave_response(chain, pump_frequency=1e9, pump_power=7e-14,
                            signal_frequency=.9e9, signal_power=7e-20)
print(10*np.log10(r.gain), r.idler_conversion)
```

The finite-circuit model includes JJ plasma response, resistive losses,
port reflections and exact-sine Fourier projection. Set `depleted=True` for
self-consistent pump depletion, or attach `jtwpa.ResonantLoad` for resonant
phase matching. Gain is a power ratio; idler conversion is a photon-flux ratio.
See the [4WM guide](docs/four-wave-jtwpa.md) for harmonic/sideband convergence,
normalization and limits. `jtwpa.quantum_noise` returns photon-normalized
scattering, vacuum/thermal bath noise, quadrature covariance and input-referred
added noise. Junction, ground and resonator resistances each have explicit
noise baths. See the [quantum-noise guide](docs/jtwpa-quantum-noise.md) and
[quantum-limit checks](docs/jtwpa-quantum-validation.md).

## KIT and phenomenological 3WM circuits

```python
import numpy as np
from twpa import Cell, Section, Design, Bias, Nonlinearity, linear, small_signal

unloaded = Cell(inductance=54e-12, capacitance=11e-15,
                stub_inductance=1.3e-9, participation=.03)
loaded = Cell(inductance=54e-12, capacitance=4e-15,
              stub_inductance=.5e-9, participation=.02)
design = Design((Section(unloaded, 23), Section(loaded, 7),
                 Section(unloaded, 24)), repeats=900,
                tail=(Section(unloaded, 5),), name="custom")
bias = Bias(pump_current=15e-6,
            nonlinearity=Nonlinearity(three_wave=380., kerr=1e6))
grid = np.linspace(0, 11e9, 32769)
response = linear(design, grid, bias.impedance)
gain = small_signal(design, bias, np.linspace(2e9, 6e9, 101), frequency=grid)
```

Use SI units. Cell capacitance is per finger, with two shunt fingers per cell.
Inductances are electrical values at the operating bias. TWPA current amplitudes
are peak, not RMS; propagation uses cell index, rad/cell and nepers/cell.
`Nonlinearity` supplies phenomenological 3WM/Kerr coefficients; omitting it
selects the KIT current law. Supplying different coefficients does not add
Josephson plasma response or full 4WM to these APIs; use `twpa.jtwpa` above.

`linear` returns Bloch dispersion, attenuation, stopbands and finite S11/S21.
`network_scattering` handles arbitrary row/turn blocks. `nonlinear` includes
pump depletion and optional spatial film gradients; `compression` sweeps input
power and reports the first bracketed 1-dB compression crossing.

## TLS loss and amplifier noise

```python
from twpa import preset, IDCGeometry, oxide_design, TLS, coupled

design, bias = preset('v2-tls')
coated, bare, fields = oxide_design(
    design, IDCGeometry(mesh=5e-9), oxide_thickness=30e-9)
tls = TLS(loss_tangent=.001, temperature=.02,
          critical_field=30., field_per_amp=2e8, linewidth=1e6)
result = coupled(coated, bias, [3.5e9, 4e9, 4.5e9], tls)
print(fields.participation, result.gain_db, result.added_noise)
```

The pipeline calculates IDC edge fields, scales circuit capacitance from
coated/bare electrostatic energy, and couples lossy propagation to amplifier
transfer and thermal bath covariance. Material and device values above are
illustrative. Verify mesh/domain convergence and calibrate the field/current
conversion for a physical device.

The quantum-noise propagation model is small signal: pump attenuation is
included, while conversion-induced depletion is handled by the separate
classical nonlinear solver. Gain assumes forward propagation; finite port
reflection feedback is not included. The JoSIM comparisons currently cover
the resonant JPA and separate JTWPA models, not this KIT/TLS loss pipeline.

## JoSIM circuit validation

JoSIM is an optional external simulator. The [backend guide](docs/josim-backend.md)
covers native and Windows/WSL execution, generated junction/SQUID circuits,
coherent port-wave measurements and reproducible run artifacts.

```powershell
python tools/validate_jpa_josim.py --executable /path/to/josim-cli `
  --wsl Ubuntu --output validation/runs/example --step 5e-13 --compression
```

Omit `--wsl` for a native executable. Tested circuit/model gain and conversion
differences are below 0.08 dB for the documented benchmark points. See the
[validation report](docs/josim-validation.md) for parameters, timestep/window
checks and limits. Classical JoSIM runs do not validate quantum added noise
or vacuum squeezing noise.

The [JTWPA circuit report](docs/jtwpa-josim-validation.md) covers uniform
and resonator-loaded 4WM chains, passive dispersion and strong-probe depletion.
Tested signal/idler differences are below 0.01 dB, with timestep, harmonic
and acquisition-window sensitivity recorded separately. Reproduce runs with
`tools/validate_jtwpa_josim.py`.

## Examples, presets and validation

| Example | Purpose |
| --- | --- |
| `resonant_jpa.py` | Kerr JPA gain/noise and effective flux-pump squeezing |
| `four_wave_jtwpa.py` | Uniform and resonator-loaded 4WM Josephson chains |
| `quantum_jtwpa.py` | Quantum-limit and internal thermal noise in 4WM chains |
| `compare.py` | Compare electrical TWPA parameter sets |
| `oxide_pipeline.py` | Electrostatics -> lossy circuit -> amplifier noise |
| `depletion_and_gradients.py` | Pump-depleted propagation with spatial profiles |

TWPA presets are example electrical designs, not universal material models.
`v1-target`, `v1-cpw`, `v1-idc` and `v2` provide alternative circuit estimates;
`v2-tls`/`v2-standard`, `v2-folded` and `v2-four-pass` use a calibrated circuit
with different device lengths. Their names encode circuit and length, not an
automatically simulated routing geometry. Configure gradient or row/turn
circuits explicitly. Preset coefficient conventions and migration details
are documented in the [source audit](docs/source-audit.md).

Tests cover analytic circuit limits, energy conservation, source regression,
loss/noise limits, JPA roots and stability, commutator preservation, squeezing,
and transient/small-signal agreement. Frequency grids must resolve all relevant
bands below circuit poles; verify numerical convergence for a new design.
Passing tests does not establish agreement with a fabricated device.

See the [Josephson model roadmap](docs/josephson-roadmap.md) for JTWPA and further circuit
implementation plans. Contributions are licensed under the [MIT License](LICENSE).
