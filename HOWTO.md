# HOWTO: simulate the main amplifier families

This is the **electrical-simulation stage** of the
[overall design cycle](docs/design-workflow.md). Use
[Tools and setup](docs/tools-and-setup.md) to choose your environment first.
The numbered examples below use the Python core; actual JoSIM runs, layout
checks and field extraction add the external tools described there. After an
electrical candidate is useful, follow the cycle through layout, checking,
extraction and revised simulation rather than assuming that a successful
electrical example is an accepted physical design.

These examples simulate superconducting microwave amplifiers. A strong
microwave pump transfers energy into a weaker signal through a nonlinear
inductance; an additional frequency, the idler, participates in the process.
The scripts calculate signal growth and, where supported, its accompanying
noise. Start with [What is being simulated?](docs/physical-overview.md) for
the physical picture and the differences between device families.

The examples are small learning cases with illustrative electrical values.
Several deliberately have modest gain: they demonstrate a calculation, not
an optimized high-gain amplifier or a fabrication-ready design. Each section
below explains the circuit, the purpose of its parameters, and how to read
the printed results. The [examples index](examples/README.md) also covers the
longer exploratory scripts outside this sequence.

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

## Reading the output

Electrical parameters use SI units: amperes, farads, henries, ohms, watts and
kelvin. Frequency inputs are in hertz, so `1e9` is 1 GHz and `0.02` kelvin is
20 mK. Effective resonator decay and nonlinear coupling rates are angular
rates in rad/s; a linewidth in Hz needs the factor `2*pi`.

Power gain `G` is output signal power divided by incident signal power, and
`gain_dB=10*log10(G)`. Unity gain is 0 dB. Josephson `idler_conversion` is
instead a photon-flux ratio: divide power by photon energy `h*f`, with
Planck constant `h` and physical frequency `f`, to count photons per second.
Multiply the idler photon ratio by the idler/signal frequency ratio for its
power ratio.

Input-referred added noise divides output fluctuations by `G` and subtracts
the original signal-input noise. The symmetrized convention assigns vacuum
one half quantum. For a phase-preserving amplifier with `G>1`, the ideal
lower bound is `(1-1/G)/2`, approaching one half only at large gain.
At coincident signal/idler frequencies, use the quadrature gains and noise
explained in example 04. See [Caves](docs/references.md#caves1982),
[Caves et al.](docs/references.md#caves2012), and the
[physical output guide](docs/physical-overview.md#how-to-read-the-results).

## KIT and TLS loss

A kinetic-inductance traveling-wave amplifier (KIT) uses a superconducting
film whose inductance changes with current. Pump and signal propagate through
many repeated circuit cells. The desired amplified signal leaves the far end,
so the reported gain is a transmission measurement. Example 01 asks how its
weak-signal prediction changes when a finite signal is allowed to take energy
from the pump. Example 02 asks how surface-dielectric loss changes gain and noise.

```powershell
python examples/howto/01_kit.py
python examples/howto/02_tls_loss.py
```

Example 01 uses `Cell`, `Section` and `Design` to build a 50-ohm
low-frequency ladder. `Cell.capacitance` is **per finger**, with two fingers
per cell. `Bias` sets peak pump current, DC bias and the kinetic-inductance
scale. This uniform line has no loading stopband for automatic pump selection,
so it supplies `pump_frequency=7e9`. Here `f_p`, `f_s` and `f_i` mean pump,
signal and idler frequencies in hertz. Three-wave mixing (3WM) obeys
`f_i=f_p-f_s`, giving a 4.2 GHz idler for the 2.8 GHz signal in example 01.
A passband is a frequency range that supports propagation in the periodic
cell; all three waves need a suitable propagating response.

`small_signal` gives undepleted forward gain. `nonlinear` propagates three
coherent currents together, including classical depletion; it does not return
quantum noise. Try changing `repeats`, `pump_current` or the finite probe
amplitude one at a time. Multiple `Section` objects make a loading pattern.

Read `weak_gain_dB` and `finite_gain_dB` as signal power gain in decibels;
0 dB is unity gain and 10 dB is tenfold power. `output_pump_peak_A` is the
remaining pump's peak current, in amperes. The passive power-balance error
measures departure from reflected plus transmitted power equaling incident
power in this lossless circuit and should be close to zero.

Example 02 connects the complete loss pipeline:

Its interdigitated capacitor (IDC) consists of interleaved metal fingers.
Strong electric fields near the finger edges store some energy in a surface
oxide. Two-level systems (TLS) in that oxide can absorb energy. The simulation
first estimates that electric-energy fraction, called participation, and then
uses it in the propagating amplifier's loss and noise calculation.

1. `oxide_design` computes coated/bare field energies and oxide participation.
2. The energy ratio rescales calibrated electrical capacitance.
3. `coupled` evolves pump attenuation, local signal/idler transfer and bath noise.

Check mesh and domain convergence and calibrate the loss tangent and
field/current conversion for your geometry. The saturation-hole linewidth
is phenomenological. This forward small-signal noise model includes pump
attenuation, but conversion does not deplete that pump and reflected waves
are not fed back into its coupled-mode equations, which evolve traveling-wave
amplitudes along the line. `participation` refers to the modeled oxide edge
domain, not a full chip. `capacitance_scale` is the coated/bare electric-energy
ratio. Compare `added_quanta` with the noise convention explained below and
inspect how the output pump decreases with loss. See [KIT/TLS](docs/kit-twpa.md) and its
[references](docs/references.md#kit-and-dielectric-loss).

## Kerr JPA

A resonant Josephson parametric amplifier (JPA) stores microwave energy in
a nonlinear resonator. Here the amplified signal is reflected from the same
port through which it entered. The calculation uses a single effective mode:
it starts with resonance and damping, rather than individual junctions or mask
dimensions. A Kerr nonlinearity means that the resonance frequency shifts
with the number of photons stored in that mode.

```powershell
python examples/howto/03_kerr_jpa.py
```

`Resonator` takes effective resonance, external/internal linewidths, Kerr
coefficient and internal temperature. Convert linewidths in Hz to angular
rates using `2*pi`; do not pass MHz directly to `kappa_external`.

`kappa_external` is the energy-decay rate into the measurement port, and
`kappa_internal` is the rate into internal losses. Their sum is the total
decay rate; the field amplitude decays at half that rate. `kerr` is the
angular resonance shift per intracavity photon. Its negative sign in the
example means that the resonance moves downward as occupation grows.

`pump_points` returns every physical steady-state root with stability.
Select a stable branch using `select_pump_point`, then call `kerr_response`.
The idler is `2*f_p-f_s`. Stay away from coincident signal and pump frequencies
for this phase-preserving added-noise calculation.
Try changing pump frequency/power and inspecting root stability first.
`kerr_transient` studies a finite classical probe; vary its duration and
discard fraction. See the [JPA guide](docs/resonant-jpa.md) and
[Eichler–Wallraff reference](docs/references.md#eichler2014).

In the output, `cavity_photons` describes the pump's stored occupation,
`selected_branch` describes its linear stability, and `gain_dB` is the weak
signal's reflected power gain. More than one pump steady state can exist;
the lower branch means lower photon occupation, not lower drive frequency.
Vary pump power while inspecting these states before interpreting a gain sweep.

## Flux-pumped JPA

```powershell
python examples/howto/04_flux_jpa.py
```

The physical pump is near twice the resonator frequency. `pair_drive` is an
effective coupling rate in rad/s, not power or flux amplitude. At degeneracy
there is one signal/idler mode: use `jpa.quadratures` for two quadrature gains
and output noise spectra. Vacuum variance is 1/2. Scalar phase-preserving
`added_noise` is intentionally NaN there.

Degeneracy means that signal and idler coincide at 5 GHz in this example.
The two quadratures are cosine-like and sine-like field components relative
to the pump. One combination can be amplified while the orthogonal one is
deamplified. `quadrature_power_gains` contains their power gains;
`output_noise_variances` contains their symmetrized spectral noise. A variance
below the vacuum reference 0.5 is squeezing in this normalization. A single
phase-preserving added-noise number would obscure that phase dependence.

Remain below the effective model's instability threshold. To relate physical
coil flux to the rate, use the `SQUIDJPA` adapter and [JoSIM workflow](docs/josim-backend.md).
The conversion is not universal for arbitrary SQUID geometries.

## 4WM JTWPA

A Josephson traveling-wave amplifier (JTWPA) replaces the film's nonlinear
inductance with a chain of Josephson junctions. This example connects 64
junctions in series with capacitance from each node to ground. Each junction
also has its own capacitance and loss resistance. The finite circuit includes
reflections at physical source and load resistors.

```powershell
python examples/howto/05_four_wave_jtwpa.py
```

`jtwpa.Chain` takes physical critical current `I_c` in amperes, junction
capacitance `C_j` in farads, ground capacitance, cell count, ports
and loss resistances. `ResonantLoad` adds identical capacitively coupled LC
shunts for resonant phase matching. Harmonic balance includes finite-port
reflections and junction capacitive dynamics. Here four-wave mixing (4WM)
obeys `f_i=2*f_p-f_s`: the 1 GHz pump couples a 0.9 GHz signal to a 1.1 GHz idler.
Harmonic balance represents the circuit motion by retained sinusoidal
components and solves their amplitudes; the [method guide](docs/numerical-methods.md)
explains the numerical controls.

Use `four_wave_response` for coherent gain and `quantum_noise` for Gaussian
fluctuations with junction, ground and resonator baths. Temperatures are K,
resistances ohms, powers W. `idler_conversion` is a photon-flux ratio;
multiply by `f_i/f_s` for the idler power ratio. Photon flux is power divided
by photon energy `h*f`, where `h` is Planck's constant.

Read `gain_dB` and `idler_photon_conversion` first, then compare
`added_quanta` with `finite_gain_limit`. The latter is the quantum lower bound
for the calculated gain, not an experimental target fixed at half a quantum.
`commutator_error` is a consistency defect in the quantum transfer map and
should be small. `maximum_phase_rad` reports junction phase excursion; it is
a model-domain diagnostic, not a measured switching threshold. The short,
lossy example has modest gain and is not expected to reach a published
high-gain device's performance. For phase matching, see
[O'Brien et al.](docs/references.md#obrien2014).

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
flux bias. The name expands to radio-frequency superconducting quantum
interference device. The applied flux shifts the equilibrium phase, giving
an asymmetric potential that enables 3WM. Screening means that the loop's
own circulating current changes its total flux. Its strength is
`beta_L=L_loop*I_c/phi_bar`, with loop inductance `L_loop`, junction critical
current `I_c`, and reduced flux quantum `phi_bar=Phi_0/(2*pi)`. The flux quantum
is `Phi_0=h/(2*e)`, where `e` is the positive elementary charge. The supported
unique-equilibrium branch has `beta_L<1`.
An isolated zero-inductance dc-SQUID is not substituted for this static 3WM element.

For your circuit, create `RFSQUID(critical_current, loop_inductance,
junction_capacitance, ...)` and then `ThreeWaveChain`. Inspect
`coefficients(element)` for normalized potential derivatives `c2`, `c3` and
`c4`: `c2` sets linear curvature, `c3` gives 3WM, while
residual Kerr and harmonics remain present. The idler is fp-fs.

`linear_chain.critical_current` means Ic*c2, a **curvature scale**, not the
physical junction's switching current. The example uses it to choose an
available input-wave power. Actual junction currents/phases come from the
finite-circuit solution. `dc_phase` reports the pump's static rectification.
`flux_Phi0` is applied flux divided by the flux quantum.
`maximum_DC_shift_rad` is the pump-induced shift of the average phase, in
radians: an asymmetric potential can rectify an oscillating drive.
See the [element potentials and limits](docs/three-wave-josephson.md) and
[Zorin's physical circuit description](docs/references.md#zorin2016).

## SNAIL 3WM

A superconducting nonlinear asymmetric inductive element (SNAIL) has a loop
with one small junction on one branch and several larger series junctions
on the other. Its asymmetry and magnetic bias enable 3WM. The example uses
three large junctions, a small/large critical-current ratio of 0.3 and a bias
of 0.35 flux quanta. A 1 GHz pump couples the 0.45 GHz signal to a 0.55 GHz idler.

```powershell
python examples/howto/07_snail.py
```

`SNAIL` takes large-junction Ic, small/large ratio alpha, number of large
series junctions, large-junction C/R and flux bias. The supported model needs
n>=2, 0<alpha<1/n and an unambiguous equilibrium. Small-junction capacitance
scales as alpha and resistance as 1/alpha. It assumes identical large JJs
and no intermediate-island ground capacitance.

Here `n` counts the large junctions and `alpha` is the small/large critical-
current ratio. `C` and `R` mean junction capacitance in farads and loss
resistance in ohms. `equilibrium_small_JJ_phase_rad` is the small junction's
static phase; the potential derivatives have the meanings given in the
rf-SQUID example. The model retains residual higher-order nonlinearity, so
a small quartic coefficient alone does not establish Kerr-free operation.
See [Frattini et al.](docs/references.md#frattini2017) for the element and
[Roudsari et al.](docs/references.md#roudsari2023) for traveling-wave context.

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

A netlist is a text description of circuit components and their connections.
JoSIM evolves the specified circuit in time, providing a classical comparison
with the Python frequency-domain models. Its Josephson elements include the
sine current law, capacitance and resistance: the resistively and capacitively
shunted junction (RCSJ) model. Flux-biased examples contain physical inductors
and coils. See the [JoSIM reference](docs/references.md#josim2019).

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
WSL means Windows Subsystem for Linux.
Use `--corrector 4` only if the executable advertises that feature. This
weak classical probe does not by itself measure gain or quantum noise.
The [backend guide](docs/josim-backend.md) and validation reports give
family-specific `tools/validate_*_josim.py` commands for pump/probe subtraction
and port-wave extraction. Pin the simulator binary throughout a study.

## Check a new design

GDS is the layout exchange file, and DRC means design-rule checking of drawn
widths, spacings and related geometric constraints. The companion manifest
records the design and process assumptions. Connectivity/DRC checks can detect
geometric mistakes; electrical gain still requires a circuit or field model.

For mask geometry, follow the [layout HOWTO](docs/layout.md): select a process
JSON, set device geometry in an [example design](examples/layout/kit_small.json),
generate GDS, then run independent geometry/connectivity and DRC checks. The
[commented all-family script](examples/layout/generate_all.py) covers KIT,
JPA, SQUID JPA, 4WM, rf-SQUID and SNAIL layouts. Electrical sizing for a bare
JPA or uniform 4WM chain rejects incompatible junction current/capacitance
targets; loop inductance and parasitics still require extraction.

For single-junction JPAs, the [London-network workflow](docs/magnetic-extraction.md)
combines openEMS capacitor extraction with FastHenry interconnect, return,
contact, kinetic energy and signed coupling extraction. The
[commented example](examples/em/london_jpa.py) prepares either a physical pickup
loop or the original galvanic-port geometry without launching a solver.
The full driver checks field convergence and the annotated full-sine circuit
against JoSIM. London L already contains kinetic inductance; add it only once.
Electric stray capacitance outside the MIM coupon and packaging remain outside
this hybrid model. See the [numerical report](docs/magnetic-jpa-validation.md).

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

The [numerical-method guide](docs/numerical-methods.md) explains how to check
numerical convergence separately from physical stability.
