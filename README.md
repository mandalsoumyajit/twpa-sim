# TWPA Sim

Parameterized traveling-wave parametric amplifier circuit and envelope solvers.
Repository: **twpa-sim**. Python import: **twpa**.

Circuit values, loading pattern, device length, port impedance, operating
currents, nonlinear coefficients, dielectric participation, TLS properties and
numerical grids are explicit inputs. The recovered NY CREATES solvers provide
the initial kinetic-inductance presets and regression references.

## Install and test

Python 3.10+, NumPy and SciPy are required. MATLAB and KLayout are not needed.

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -e .
.\.venv\Scripts\python -m unittest discover -s tests -v
.\.venv\Scripts\python examples\compare.py
.\.venv\Scripts\python examples\oxide_pipeline.py
.\.venv\Scripts\python examples\depletion_and_gradients.py
```

## Define a design

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
# Extracted current-envelope coefficients, independent of a KIT current law:
bias = Bias(pump_current=15e-6,
            nonlinearity=Nonlinearity(three_wave=380., kerr=1e6))
grid = np.linspace(0, 11e9, 32769)
response = linear(design, grid, bias.impedance)
gain = small_signal(design, bias, np.linspace(2e9, 6e9, 101), frequency=grid)
```

All inputs use SI. Capacitance is **per finger**, with two shunt fingers per
cell. Inductances are electrical values at the operating bias; they do not
automatically rescale when bias changes. Currents are peak, not RMS. Wavenumber
and amplitude attenuation are rad/cell and nepers/cell; position is cell index.

`Nonlinearity` supplies phenomenological 3WM and Kerr coefficients in A^-1 and
A^-2. If omitted, `Bias` derives KIT coefficients using the scaling and DC
currents. This interface supports devices that fit the ladder and three-mode
envelope approximation. Josephson plasma dynamics and full 4WM sideband
generation require additional microscopic models; generic naming does not
imply those models are already implemented.

Sections may have odd or zero counts. A period repeats and appends optional
tail sections. Cell topology accepts `shunt-series` or `symmetric` (the
presentation's centered circuit). `network_scattering` accepts arbitrary
repeated blocks for distinct row/turn assemblies.

## Solvers

| API | Purpose |
| --- | --- |
| `linear` | Bloch phase, attenuation, stopbands and finite S11/S21 |
| `network_scattering` | Nonperiodic row/turn/block scattering |
| `small_signal` | Analytic undepleted-pump 3WM gain and pump search |
| `constant_langevin` | Original constant-pump lossy amplifier gain and bath noise |
| `coupled` | Distributed TLS-dependent pump attenuation, transfer and bath covariance |
| `nonlinear` | Full complex three-current CMEs with depletion and Kerr phases |
| `compression` | Power sweep and first bracketed 1-dB compression crossing |
| `route_profile` | Configurable linear, sine, folded and spiral WKB coordinates |
| `idc_fields` | Recovered finite-volume IDC edge-field solver |
| `oxide_design` | Computed coated/bare energy ratio and participation -> circuit |
| `electrostatic` | General rectangular dielectric/electrode-map solver |

Gradient phase is integrated as an extra nonlinear ODE state. Seven sine cycles
model the original meander x-coordinate; one/two folded cycles model two/four
passes. Inductance sensitivity uses a configurable centered difference.
Nonlinear coefficients remain spatially uniform, as in the original scripts.

## Coupled TLS loss + amplifier

```python
from twpa import preset, IDCGeometry, oxide_design, TLS, coupled

design, bias = preset('v2-tls')
coated, bare, fields = oxide_design(
    design, IDCGeometry(mesh=5e-9), oxide_thickness=30e-9)
# Illustrative values; replace with process characterization.
tls = TLS(loss_tangent=.001, temperature=.02,
          critical_field=30., field_per_amp=2e8, linewidth=1e6)
result = coupled(coated, bias, [3.5e9, 4e9, 4.5e9], tls)
print(fields.participation, result.gain_db, result.added_noise)
```

This connects the source's Laplace geometry and energy-ratio prescription to
lossy ABCD dispersion and Langevin amplifier equations. Participation and
capacitance ratios are calculated, replacing the old plot's hardcoded tables.
`IDCGeometry` controls finger width, gap, metal thickness, permittivities,
domain and mesh. The default 2.5 nm mesh reproduces the original boundaries.
Demonstrate mesh/domain convergence for new geometries. Half-gap domain energy
ratios scale calibrated circuit capacitances. Unlike cross-sections need
separate extraction. The helper replaces cell participation with the coating
participation and does not add an independent residual surface-loss channel.

TLS loss uses `C*(1-i*p*tan_delta)` with `exp(+i*omega*t)`, thermal
`tanh(h*f/(2*kB*T))`, and field saturation. `TLS.from_density` accepts spectral
density in J^-1 m^-3 and dipole in C m. Local peak `field_per_amp` requires
calibration. A phenomenological Lorentzian saturation hole is an explicit
extension, not a recovered microscopic driven-TLS calculation.

`coupled` evolves pump attenuation, accumulated mismatch, photon-normalized
signal/idler transfer and thermal bath covariance. Added noise is input-referred
symmetrized quanta: signal-input noise is excluded, input-idler vacuum and
internal thermal baths are included. Lossless added noise is `(1-1/G)/2`;
a zero-temperature passive attenuator gives `(1-eta)/(2*eta)`.

Original-source compatibility uses `v2-tls`, `pump_loss=False`, constant effective
loss at zero temperature and the source's pump frequency. `constant_langevin`
directly reproduces the rotating-frame solution with a matrix-exponential
covariance integral replacing the original trapezoidal quadrature. Regression
fixtures cover all six source loss cases.

The quantum-noise solver is small signal: pump attenuation is included but
conversion-induced depletion is not. `nonlinear` separately models depletion
and compression without quantum noise. Envelope gain uses forward propagation
and homogenized periodic coefficients; reflection feedback and distinct tail
nonlinearities are excluded. Tail cells enter gain length and exact scattering.

## Presets and provenance

| Preset | Source and total cells |
| --- | --- |
| `v1-target`, `v1-cpw`, `v1-idc` | MATLAB electrical models, 49,361 |
| `v2` | Historical MATLAB model, 56,700; not the final layout |
| `v2-tls`, `v2-standard` | Recovered calibrated TLS circuit, 49,950 |
| `v2-folded` | Same circuit, report's two-column length, 41,472 |
| `v2-four-pass` | Same circuit, report's four-pass length, 39,204 |

The TLS source uses `eps2=2/(I_star^2+I_dc^2)`; other sources use `1/(...)`.
`Bias.kerr_factor` retains the distinction, with `2` for TLS/layout presets.
The layout names encode circuit and length only; supply routing profiles or
row/turn circuits explicitly.

See [source audit](docs/source-audit.md). Source hashes and captured results
are in `tests/fixtures/original_solvers.json`. Regenerate using:

```powershell
python tools/capture_reference.py PATH\TO\NY_CREATES\scripts
```

Original files remain intact. GDS, correspondence, plotting side effects and
machine-specific output paths are not required. Regression and analytic tests
check numerical behavior, not measured-device accuracy.

## Numerical behavior and publication

Frequency grids start at DC, increase strictly and must resolve each band below
effective-stub poles. Phase uses the full supercell trace and continuation from
DC, replacing whole-device transmission phase and fixed 8.46 GHz branching.
Scaled matrix powers retain S21 in dB through long stopbands. Check grid
convergence; a phase-step guard cannot detect an entirely skipped band.

Pump search spacing, offsets and band index are configurable. Supply the pump
explicitly without a stopband. Analytic gain marks unsupported frequencies NaN;
ODEs reject them. Compression returns NaN for an unbracketed crossing.

The repository is hosted privately at
[mandalsoumyajit/twpa-sim](https://github.com/mandalsoumyajit/twpa-sim), with
packaging, examples and GitHub CI configuration. A redistribution license is
not assumed; choose an appropriate license before a public release.

The [Josephson model scope](docs/josephson-roadmap.md) covers resonant JPA,
4WM/3WM JTWPA models, and an optional JoSIM transient-simulation backend.
Those extensions are planned, not yet implemented.
