# Layout back-annotation across amplifier families

The purpose of back-annotation is to replace electrical assumptions with
parameters tied to a particular checked mask. A schematic can connect two
components with a perfect wire; the drawn connection stores magnetic and
kinetic energy. In a SQUID or SNAIL, that energy also changes flux screening
and the junction phases. In a traveling-wave circuit, different end
capacitors and finite return paths can change reflection and conversion.

Use this page after the [electrical HOWTO](../HOWTO.md) and
[layout checks](layout.md). The organizing loop remains **simulate → generate
layout → check → extract → back-annotate → compare → revise**. Read
[Tools and setup](tools-and-setup.md) before launching external programs.

## Two different routes

| Layout | Returned circuit | Physical parameters supplied by the route |
| --- | --- | --- |
| SQUID JPA | `LayoutNetwork` | Each junction, its actual area, separate shunt capacitor, full signed metal inductance matrix and physical terminals |
| 4WM JTWPA | `LayoutNetwork` | Every junction and shunt/end capacitor; finite signal/return conductors and all matrix couplings |
| rf-SQUID / SNAIL 3WM JTWPA | `LayoutNetwork` | The same finite network, including loop arms and every large SNAIL junction; applied flux is declared separately |
| KIT | `KITPeriodicModel` | Drawn cell counts, sheet-squares inductance and turn matching, plus explicitly supplied capacitance/stub/material calibration |

The dedicated [single-JJ JPA adapter](magnetic-extraction.md) remains available
for galvanic or pickup-loop circuits. Its published field-convergence and
JoSIM benchmark is separate from the new general-network evidence. See the
[new validation report](layout-network-validation.md) for accepted and
rejected runs. An implemented family adapter is not a certificate that an
arbitrary mesh, material profile or amplifier design has been validated.

The KIT route is a **calibration-backed periodic reduction**. It does not
extract the entire serpentine with openEMS or FastHenry. Its record retains
that distinction explicitly; see the KIT section below.

## Josephson route: keep the finite network

1. Choose `JosephsonLayout(family=...)`, generate GDS and keep its manifest.
   The supported family names are `squid_jpa`, `jtwpa4wm`, `rfsquid3wm` and
   `snail3wm`. Run geometry/connectivity/DRC checks before extraction.
2. Call `layout_magnetic_job` with an explicit London stack and fixed terminal
   patches. Run it with the **Whiteley superconducting FastHenry** adapter.
3. Extract each `CAP` component using `capacitor_job(..., component_id=...)`
   and **openEMS/CSXCAD**. Refine mesh, domain and duration before certification.
4. Refine the magnetic bulk grid, thickness discretization and local contact
   grid, and compare the fixed current basis at high frequency. Call
   `certify_magnetic`; a failed result is retained and rejected for annotation.
5. Pass the complete magnetic certificate and a capacitor-ID-to-certificate
   dictionary to `back_annotate_layout`. Simulate its returned `LayoutNetwork`.
6. Generate the same physical circuit with `circuit.netlist(...)`, run
   **JoSIM**, and compare passive scattering, signal, idler and junction phases.
   Refine circuit harmonics, sidebands, timestep and acquisition independently
   of the field mesh. Then compare the response with the original design goals.

The [commented preparation example](../examples/em/layout_network.py) launches
no external solvers. It lets you inspect the terminal map and capacitor coupons
first. For a native extraction/comparison run, use a fresh output directory:

```powershell
python tools/validate_layout_network.py --family rfsquid3wm `
  --output validation/runs/rf-iteration-001 `
  --installation fasthenry-installation.json `
  --openems-python C:/path/to/openems/python.exe `
  --openems-dll-directory C:/path/to/openEMS `
  --executable /path/to/josim-cli --wsl Ubuntu --corrector 4
```

Replace paths with your installed tools; omit `--wsl` for native JoSIM.
Corrector 4 requires a JoSIM build advertising that capability. The driver
uses float64 NPY traces by default, which also require a compatible build.
`--output-format csv` supports upstream CSV output, but timestamp/voltage
precision must pass the same numerical comparisons. `--cells`
sets the chain length; `--large-junctions` sets the SNAIL series-arm count.
Supply `--gds` to use a previously generated GDS/manifest pair. The default
one-cell cases exercise finite circuit assembly and conversion, not useful
traveling-wave gain or bandwidth.

The electric driver samples 1–100 GHz, including the 5 GHz pump and the
2/3 GHz or 4/6 GHz signal/idler pairs used in these tests. Each capacitor
record retains its extraction band. A new operating band needs an appropriate
fit and checks; full-chip electric coupling is still omitted.

The driver defaults to one FastHenry worker, an 8 GiB address-space cap,
900 s per native job, a 2 um bulk pitch, 0.125 um local pitch, two thickness
subdivisions and 1 um terminal patches. These are starting controls, not a
promise of convergence. `--pitch-um`, `--local-pitch-um`, `--nh`,
`--patch-side-um`, `--timeout` and `--memory-gib` expose those choices.
`--contact-control-pitch-um` adds an intermediate coarser contact grid after
retaining earlier rejected levels. This changes the declared comparison,
not the numerical acceptance limits. The report lists every grid used.
`--dense-threads` bounds dense-library threads to 1–4.
The current film adapter requires axis-aligned (Manhattan) two-metal geometry.
Its Cartesian basis rejects jobs above 40,000 filaments.
Large chains/SNAIL layouts may need a better mesher before full extraction is
practical. Increasing a timeout does not fix a mesh-convergence failure.

A patch is a fixed finite electrode at a physical reference plane. Changing
its size changes the physical partition between the field network and the
lumped component. Keep it fixed during convergence comparisons. Refine the
contact grid rather than relabeling different patches as the same problem.

The compact SNAIL benchmark uses two large junctions, a 40 um cell pitch,
a 16 um branch height, 8 by 8 um full-capacitor dimensions, 1 um wires and
0.5 um contacts. Prepare it with
`examples/em/layout_network.py --family snail3wm --compact-snail --output ...`.
For its native driver, supply that GDS/manifest with `--gds` and use
`--patch-side-um .5 --local-pitch-um .125 --contact-control-pitch-um .17`.
The rf-SQUID benchmark uses `.17`/`.25` local grids with 1 um patches; the
SQUID/4WM benchmarks use `.125`/`.17` with 1 um patches. Other dimensions
and mesh levels require their own convergence study.

The earlier 2 um-wire/1 um-contact SNAIL failed its fixed-terminal coupling
check. The narrower design was regenerated, rechecked and fully re-extracted;
its inductance was not scaled from the rejected layout. Its identical local
capacitor coupons may reuse the native electric evidence after exact geometry
comparison.

The compact SNAIL electric test uses `--electric-port-ohm 10`,
`--electric-port-control-ohm 20` and `--electric-timesteps 200000`. These
are numerical measurement-port/settling controls, not added device losses.
The 1 ohm low-frequency case and the shorter transient are retained as
rejected numerical attempts. The validation report gives the measured
capacitance sensitivity to these controls.

### What is in the circuit?

M1 and M2 name the two metal layers, and VIA denotes an actual filled contact.
Each junction and capacitor has a terminal patch on both layers. Input/output
ports use separate signal and return patches; the return is a finite drawn
conductor. A spanning star connects the patches within each physical conductor
as a **current-coordinate basis**. Its complete signed inductance matrix is
retained, so the star does not merge the conductor into an ideal circuit node.
Negative mutual entries keep their native port orientation.

The circuit has one full sine-law RCSJ element per junction. RCSJ means
resistively and capacitively shunted junction. Junction critical current
`I_c` and capacitance `C_j` equal the measured drawn window area times the
explicit process densities. Each separate openEMS shunt capacitance `C_s` is
added once, independently of `C_j`. Values are in amperes, farads and ohms.
`IND` annotations do not add another inductor: their drawn wire is already in
the extracted matrix. Total London inductance includes superconducting kinetic
energy; adding another sheet-kinetic inductance would double count it.

The general four-control magnetic certificate covers **total inductance**.
It does not claim the separate kinetic-energy derivative and PEC-limit
energy decomposition reported by the dedicated JPA benchmark.

The capacitor model is a cropped PEC component extraction. The metal network
is a London magnetoquasistatic extraction. Their combination omits full-chip
electric cross-capacitance, radiation and packaging. Constant junction and
dielectric loss resistances are explicit circuit assumptions; they are not
measurements derived from PEC or lossless London calculations.

The nodal equation makes the partition explicit:

```text
C * d²phi/dt² + G * dphi/dt
  + B_L * inverse(L) * (transpose(B_L) * phi - f)
  + B_J * [I_c * sin(transpose(B_J) * phi / PHIBAR)] = I_drive(t)
```

Here `phi` is the vector of node fluxes in webers; its time derivative is
node voltage. `C` is the assembled capacitance matrix in farads, `G` is
conductance in siemens (including terminations), and `L` is the full signed
metal-inductance matrix in henries. `B_L` and `B_J` are oriented connection
matrices for metal branches and junctions: a column has +1 at its positive
terminal and −1 at its negative terminal. `f` is the imposed branch-flux
vector in webers, `I_c` is the vector of junction critical currents in amperes,
and `I_drive` is injected current. `PHIBAR = Phi0/(2*pi)` converts node flux
to phase; the sine acts separately on each junction. The Python implementation
uses `theta = phi/PHIBAR` as its dimensionless coordinate. JoSIM solves the
same RCSJ/inductor/capacitor network in time.

### Flux and nonlinear simulation

`circuit.cycles` lists fundamental superconducting loops, using the ordering
of metal branches followed by junctions. Capacitors and resistors do not close
a superconducting DC loop. Each row has the closing junction oriented
positively. `flux_bias` is one signed applied flux per loop in units of the
flux quantum, `Phi0 = h/(2e)`. A scalar applies the same value to each loop.
Inspect this basis before assigning a spatially varying bias.

`equilibrium(flux_bias)` solves static junction phases with geometric
screening. Its initial guess selects a local branch; a positive static
Hessian does not establish global uniqueness or simulate switching between
metastable states. `periodic(pump_frequency, pump_power, ...)` solves the DC shift and
pump harmonics together. Frequency is in hertz; power is available incident
power in watts at the input impedance. `flux_pump` adds a sinusoidal applied
flux, also in units of `Phi0`. The JoSIM deck realizes the same prescribed
flux using explicit coupled ideal current-driven coils. **No bias-coil mask
has been extracted.** The coil normalization is arbitrary; its inductance
matrix is chosen positive definite and its mutuals produce the declared flux.

`passive(frequencies, flux_bias=...)` returns complex scattering coefficients
at every measurement port, for excitation at input port zero by default.
`response(point, frequencies, sidebands=...)` linearizes the full junction
sines about the periodic pump. Its `waves` array is indexed by signal sample,
sideband, then output port. Sideband `k` has signed frequency `f_s + k*f_p`.
The idler uses `k=-2` for 4WM or `k=-1` for 3WM and conjugation of the
negative-frequency result. Convert power-wave amplitudes to photon amplitudes
with `sqrt(f_s/f_i)` for the idler. The validation driver performs that
conversion and fits actual input/output voltages from JoSIM traces.

A small residual and harmonic/sideband refinement establish numerical
self-consistency. The DC Hessian and AC phase guard do **not** certify
Floquet/dynamical pump stability. The general extracted-network response is
classical; the existing uniform-chain quantum APIs must not be applied to it
by silently replacing its nonuniform parameters with averages. See
[numerical methods](numerical-methods.md) and
[quantum-noise conventions](jtwpa-quantum-noise.md).

### Reuse and revision

`reuse_capacitor(gds, component_id, certified)` reuses a native electric result
only for identical translated electrode polygons, stack and numerical settings.
It stores the complete sealed parent evidence and binds the new component to
its own GDS/manifest hashes. Distinct end geometry requires a new extraction;
a half-height capacitor is not assigned half the full capacitor's value.
The CLI accepts `--electric-library` certificates and `--magnetic-primary`
or `--magnetic` for exact compatible reuse. Geometry, process, design,
manifest or terminal changes invalidate the corresponding magnetic evidence.

If a convergence gate fails, refine that numerical control and retain the
failed record. If converged Python/JoSIM results agree but the design misses
gain, bandwidth or loading targets, revise geometry or drive, regenerate,
recheck and repeat the affected extractions. Do not fit the extracted values
to force agreement with a preferred circuit curve.

## KIT route: geometry plus explicit calibration

The [commented KIT example](../examples/em/kit_backannotation.py) generates a
small checked serpentine, binds a calibration and returns the existing
propagation model. It needs Python and KLayout, with no external field program.
Its demonstration capacitance and stub values come from the documented
research reference; changing a finger length does not make those values valid
for the new geometry automatically.

`kit_calibration` requires two capacitances per finger, one for unloaded and
one for loaded cells; each circuit cell contains two fingers. It also records
operating temperature, source description, zero-bias sheet inductance,
optional stub inductances, and TLS participation factors. Participation is the
fraction of electric energy assigned to the specified lossy dielectric.
A stub inductance represents the finger's series inductive correction.
Supply new measured or field-derived calibration when the geometry changes.
Temperature records the stated calibration condition; no temperature law is
inferred from the process JSON. The example's 1% participation and TLS loss
tangent are labeled demonstration assumptions, alongside the reference C/stub
values.

`back_annotate_kit(gds, calibration, bias)` uses the actual row/period/turn
counts. With sheet inductance `L_sheet`, cell pitch `p`, spine width `w` and
current ratio `r = I_dc/I_star`, the biased straight-cell inductance is
`L_cell = L_sheet*(p/w)*(1+r*r)`. `I_star` is the kinetic-inductance scaling
current; `I_dc` is applied DC current. This quadratic constitutive law is an
explicit modeling assumption. The requested electrical ratio must match the
layout's turn-design bias.

The required policy `turn_reduction='sheet_squares_matched_periodic'` checks
the original generator's tapered-width turn matching: the sum of local
biased turn inductances equals the same number of straight cells. The record
keeps the full turn-inductance distribution and the matching error. The
returned model replaces those turns by the declared periodic equivalent.
That preserves total biased inductance, but does not preserve detailed bend
dispersion, spatial nonlinear current or turn-specific capacitance.

The model exposes `linear`, `small_signal` and `coupled`. The last retains the
existing coupled TLS-loss/amplifier/noise pipeline, using the explicitly
supplied participation and loss parameters. `load_circuit` reloads a sealed
back-annotation record. Launches, package connections to separate ground
combs, nonlocal magnetic coupling and full serpentine EM remain outside this
reduction. The NY CREATES profile is a research reference, not a qualified PDK.

The physical background is covered by the [KIT guide](kit-twpa.md),
[rf-SQUID/SNAIL guide](three-wave-josephson.md) and their
[primary references](references.md). The signed finite-terminal treatment
follows the [Whiteley documentation](references.md#whiteley_fasthenry);
the electric component workflow follows
[openEMS](references.md#openems).
