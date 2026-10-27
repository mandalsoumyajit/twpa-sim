# HOWTO: simulate the main amplifier families

Start with small, commented examples in a source checkout. All device values
are illustrative. Change physical parameters deliberately and check numerical
convergence before using a result quantitatively.

## Install and run

Follow [installation](README.md#install). Below, `python` means the interpreter
where you installed `twpa-sim`; alternatively use its full `.venv` path.

```powershell
python examples/howto/03_kerr_jpa.py
```

This prints a stable pump state, signal/idler frequencies, gain and noise.
No simulator is needed. Every script has a `main()` function and comments
beside the parameters you are likely to change.

| Script | What you learn |
| --- | --- |
| [01_kit.py](examples/howto/01_kit.py) | Build a ladder; compare weak and finite probes |
| [02_tls_loss.py](examples/howto/02_tls_loss.py) | IDC participation → coupled TLS gain and noise |
| [03_kerr_jpa.py](examples/howto/03_kerr_jpa.py) | Select a stable Kerr branch; calculate reflection gain |
| [04_flux_jpa.py](examples/howto/04_flux_jpa.py) | Use quadratures at degenerate flux-pumped operation |
| [05_four_wave_jtwpa.py](examples/howto/05_four_wave_jtwpa.py) | Finite JJ chain, 4WM and quantum internal loss |
| [06_rf_squid.py](examples/howto/06_rf_squid.py) | Flux screening, DC rectification and rf-SQUID 3WM |
| [07_snail.py](examples/howto/07_snail.py) | SNAIL potential, gain and small-signal quantum noise |
| [08_josim.py](examples/howto/08_josim.py) | Generate a physical circuit; optionally run JoSIM |

## KIT and TLS loss

```powershell
python examples/howto/01_kit.py
python examples/howto/02_tls_loss.py
```

Example 01 uses `Cell`, `Section` and `Design` to build a 50-ohm
low-frequency ladder. `Cell.capacitance` is **per finger**, with two fingers
per cell. `Bias` sets peak pump current, DC bias and the kinetic-inductance
scale. This uniform line has no loading stopband for automatic pump selection,
so it supplies `pump_frequency=7e9`. Check passbands at fs, fp and fi=fp-fs.

`small_signal` gives undepleted forward gain. `nonlinear` propagates three
coherent currents together, including classical depletion; it does not return
quantum noise. Try changing `repeats`, `pump_current` or the finite probe
amplitude one at a time. Multiple `Section` objects make a loading pattern.

Example 02 connects the complete loss pipeline:

1. `oxide_design` computes coated/bare field energies and oxide participation.
2. The energy ratio rescales calibrated electrical capacitance.
3. `coupled` evolves pump attenuation, local signal/idler transfer and bath noise.

Check mesh and domain convergence and calibrate the loss tangent and
field/current conversion for your geometry. The saturation-hole linewidth
is phenomenological. This forward small-signal noise model includes pump
attenuation, but conversion does not deplete that pump and reflected waves
are not fed back into its CMEs. See [KIT/TLS](docs/kit-twpa.md) and its
[references](docs/references.md#kit-and-dielectric-loss).

## Kerr JPA

```powershell
python examples/howto/03_kerr_jpa.py
```

`Resonator` takes effective resonance, external/internal linewidths, Kerr
coefficient and internal temperature. Convert linewidths in Hz to angular
rates using `2*pi`; do not pass MHz directly to `kappa_external`.

`pump_points` returns every physical steady-state root with stability.
Select a stable branch using `select_pump_point`, then call `kerr_response`.
The idler is 2*fp-fs. Stay away from fs=fp for phase-preserving added noise.
Try changing pump frequency/power and inspecting root stability first.
`kerr_transient` studies a finite classical probe; vary its duration and
discard fraction. See the [JPA guide](docs/resonant-jpa.md) and
[Eichler–Wallraff reference](docs/references.md#eichler2014).

## Flux-pumped JPA

```powershell
python examples/howto/04_flux_jpa.py
```

The physical pump is near twice the resonator frequency. `pair_drive` is an
effective coupling rate in rad/s, not power or flux amplitude. At degeneracy
there is one signal/idler mode: use `jpa.quadratures` for two quadrature gains
and output noise spectra. Vacuum variance is 1/2. Scalar phase-preserving
`added_noise` is intentionally NaN there.

Remain below the effective model's instability threshold. To relate physical
coil flux to the rate, use the `SQUIDJPA` adapter and [JoSIM workflow](docs/josim-backend.md).
The conversion is not universal for arbitrary SQUID geometries.

## 4WM JTWPA

```powershell
python examples/howto/05_four_wave_jtwpa.py
```

`jtwpa.Chain` takes physical Ic, Cj, ground capacitance, cell count, ports
and loss resistances. `ResonantLoad` adds identical capacitively coupled LC
shunts for resonant phase matching. Harmonic balance includes finite-port
reflections and JJ plasma response. Here fi=2*fp-fs.

Use `four_wave_response` for coherent gain and `quantum_noise` for Gaussian
fluctuations with junction, ground and resonator baths. Temperatures are K,
resistances ohms, powers W. `idler_conversion` is a photon-flux ratio;
multiply by fi/fs for the idler power ratio.

Try `cells=128` and vary pump power. Check `harmonics=5, quadrature=24`
against the defaults. For a stronger coherent probe, set `depleted=True`
and increase `sidebands`. A small residual or phase guard passing does not
certify stability. See [4WM](docs/four-wave-jtwpa.md) and
[quantum noise](docs/jtwpa-quantum-noise.md).

## rf-SQUID 3WM

```powershell
python examples/howto/06_rf_squid.py
```

This is a one-junction **rf-SQUID** with a parallel loop inductance and DC
flux bias. Screening is included; beta_L<1 gives its supported unique branch.
An isolated zero-inductance dc-SQUID is not substituted for this static 3WM element.

For your circuit, create `RFSQUID(critical_current, loop_inductance,
junction_capacitance, ...)` and then `ThreeWaveChain`. Inspect
`coefficients(element)` for potential c2, c3 and c4: c3 gives 3WM, while
residual Kerr and harmonics remain present. The idler is fp-fs.

`linear_chain.critical_current` means Ic*c2, a **curvature scale**, not the
physical junction's switching current. The example uses it to choose an
available input-wave power. Actual junction currents/phases come from the
finite-circuit solution. `dc_phase` reports the pump's static rectification.
See the [element potentials and limits](docs/three-wave-josephson.md).

## SNAIL 3WM

```powershell
python examples/howto/07_snail.py
```

`SNAIL` takes large-junction Ic, small/large ratio alpha, number of large
series junctions, large-junction C/R and flux bias. The supported model needs
n>=2, 0<alpha<1/n and an unambiguous equilibrium. Small-junction capacitance
scales as alpha and resistance as 1/alpha. It assumes identical large JJs
and no intermediate-island ground capacitance.

The potential neglects geometric loop inductance. JoSIM uses a small explicit
loop with a size guard; verify that limit for a new design. Start with
`example_chain(..., dispersion_engineered=True)`, then replace the values
with extracted circuit parameters. This example is not a published device
reproduction or a universally Kerr-free element.

For rf-SQUID and SNAIL, `three_wave_response(..., depleted=True)` supports
classical depletion with additional sidebands. Quantum noise uses the
**pump-only** Floquet sector. Its first quantum frequencies are signal and
idler, unlike 4WM's mode list, which also includes pump modes. Consult
[conventions](docs/conventions.md) before comparing noise matrices.

## JoSIM circuits

Dry runs do not require a simulator:

```powershell
python examples/howto/08_josim.py --family snail
python examples/howto/08_josim.py --family rf-squid --write-netlist validation/runs/rf-example.cir
```

To execute, provide a native or Windows/WSL executable and a fresh directory:

```powershell
python examples/howto/08_josim.py --family snail --executable /path/to/josim-cli `
  --wsl Ubuntu --run-directory validation/runs/snail-probe
```

Omit `--wsl` for native execution. Saved deck filenames must also be fresh.
Use `--corrector 4` only if the executable advertises that feature. This
weak classical probe does not by itself measure gain or quantum noise.
The [backend guide](docs/josim-backend.md) and validation reports give
family-specific `tools/validate_*_josim.py` commands for pump/probe subtraction
and port-wave extraction. Pin the simulator binary throughout a study.

## Check a new design

1. Inspect passive S11/S21 and passbands. Resolve poles and loading features
   before choosing pump, signal and idler.
2. Reduce the probe amplitude and compare normalized gain in the linear regime.
3. Refine the relevant controls: KIT grids/ODE tolerances, IDC mesh/domain,
   JJ harmonics/sidebands/quadrature, or JoSIM timestep/settling/acquisition.
4. Check conservation, phase guards, roots/residuals and noise diagnostics.
   The phase-preserving bound at G>1 is (1-1/G)/2, not exactly 1/2 at finite G.
   At degeneracy use an appropriate quadrature representation.
5. Record physical parameters, branches, temperatures, controls, code version
   and simulator provenance. Compare differences with numerical sensitivity.

See [conventions and troubleshooting](docs/conventions.md),
[validation reports](docs/README.md#validation-and-provenance) and
[key references](docs/references.md). CI exercises all eight tutorial scripts;
actual JoSIM execution remains opt-in.
