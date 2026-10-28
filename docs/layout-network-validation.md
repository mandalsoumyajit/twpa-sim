# Multi-family layout-network validation

All four Josephson benchmark families passed fixed-terminal magnetic refinement,
individual capacitor extraction checks, and native JoSIM circuit comparisons.
These are specific illustrative finite circuits, with one traveling-wave cell
per chain fixture. They validate assembly, screening and conversion; they do
not demonstrate useful high gain in a long amplifier or qualify a fabrication
process. The separate KIT case checks an explicitly calibrated periodic reduction.

The [portable record](../validation/layout-network-reference.json) contains
the full signed matrices, numerical controls, component/drive data, native
self-contained netlists, tool versions, input hashes and rejected attempts.
It is a report, not a certificate that can be attached to another GDS. The
[back-annotation guide](layout-backannotation.md) explains the physical model,
symbols, external tools and the simulate → layout → check → extract → compare
→ revise workflow.

## What was extracted?

Each junction and capacitor retains two surface terminals. Finite input,
output and return conductors remain in the London magnetoquasistatic (MQS)
network. JJ denotes a Josephson junction. A spanning-star
current basis keeps every signed mutual entry; it does not merge real metal
into ideal nodes. Junction critical current and capacitance come from actual
window areas and the explicitly illustrative process densities. Internal loss
resistances are circuit assumptions, rather than fitted field-solver losses.

| Fixture | JJs / capacitors | Terminal nodes / metal branches | Native current filaments | Local grid: primary / comparison (um) |
| --- | --- | --- | --- | --- |
| SQUID JPA | 2 / 1 | 8 / 6 | 28,516 | 0.125 / 0.17 |
| 4WM JTWPA | 1 / 2 | 10 / 7 | 28,240 | 0.125 / 0.17 |
| rf-SQUID 3WM | 1 / 2 | 10 / 8 | 21,816 | 0.17 / 0.25 |
| SNAIL 3WM | 3 / 2 | 14 / 10 | 24,076 | 0.125 / 0.17 |

Every final GDS passed the independent native KLayout Ruby DRC deck with zero
violations. The SNAIL has one 1 by 1 um small junction, two 2 by 2 um large
junctions, 1 um wires and 0.5 um contacts. Its cell pitch is 40 um, branch
height 16 um and full-capacitor dimensions 8 by 8 um. The other fixtures use
the generator defaults. All material stacks and JJ process densities are
illustrative. The original KIT profile is a research reference.

## Magnetic refinement

The primary bulk pitch is 2 um and the comparison is 3 um. Film thickness
uses two subdivisions versus one. Each physical patch stays fixed: 1 um
square for SQUID/4WM/rf-SQUID, 0.5 um square for SNAIL. The native current
basis is also compared at 1 MHz and 100 GHz. This frequency check tests the
lossless London operator; it does not add retardation to the MQS model.

The limits are 2% self-inductance change, 5% relative change for strong mutuals,
and 0.005 normalized off-diagonal change. A strong mutual has
`abs(L_ij)/sqrt(L_ii*L_jj) >= 0.01`. The normalized change is
`abs(delta_L_ij)/sqrt(L_ii*L_jj)`, using the primary self inductances. These
`L_ii` and `L_jj` are branch self inductances, `L_ij` is signed mutual
inductance, and `delta_L_ij` is comparison minus primary. These tests
are sensitivities to the declared numerical controls, not rigorous bounds on
all omitted physics. Frequency dependence must remain below 1e-7.

| Fixture | Largest self change | Largest strong-mutual change | Largest normalized off-diagonal change | Frequency-matrix change |
| --- | --- | --- | --- | --- |
| SQUID JPA | 0.645% | 0.837% | 0.003579 | 2.66e-15 |
| 4WM JTWPA | 0.958% | 1.017% | 0.004944 | 3.80e-14 |
| rf-SQUID 3WM | 0.950% | 1.886% | 0.004901 | 1.77e-14 |
| SNAIL 3WM | 0.626% | 1.280% | 0.002892 | 4.82e-15 |

Total London inductance includes kinetic energy. This four-control general
certificate does not claim a separate kinetic-energy/PEC-limit decomposition.
Native reciprocity, positive matrix energy and lossless resistance gates also
passed. Changing the bounded dense backend from two to four threads changed
a fixed reference matrix by only 1.94e-15 relative; no matrices were symmetrized
or converted to unsigned couplings.

## Individual electric components

openEMS/CSXCAD sampled 1–100 GHz, including the actual 2–6 GHz pump/probe
frequencies. Coupons retain actual electrode polygons and feed necks. Only
identical translated coupons reuse native evidence; distinct end geometry is
extracted separately. The fit, tail decay and negative-conductance limits
remain unchanged. Mesh, domain and duration sensitivity must each stay below 2%.

| Actual capacitor marker | Extracted C (fF) | Largest refinement change |
| --- | --- | --- |
| 12 by 12 um shunt | 178.534059 | 0.2029% |
| 12 by 6 um end | 97.321083 | 0.3732% |
| 8 by 4 um end | 51.620205 | 1.7123% |

`C` denotes capacitance; one femtofarad (fF) is 1e-15 farad.
The half-height 12 by 6 um capacitor is not assigned half the full shunt
value: fringing and clipped feed geometry are retained. Junction capacitance
is added separately, once per junction.

For the smaller SNAIL coupon, the original 1 ohm extraction port failed the
low-frequency numerical-conductance test. A 10 ohm measurement port and
200,000 timesteps passed, with independent comparisons at 20 ohms and
100,000 timesteps. These source/settling controls are unrelated to the
10,000 ohm device loss resistances used in the circuit. The measured source-
impedance sensitivity was 0.00139%; the shorter 50,000-step attempt remains rejected.

## Native circuit comparison

Python solves the full junction sines, DC rectification, pump harmonics and
Floquet sidebands (frequencies shifted by integer pump frequencies). RCSJ
means a resistively and capacitively shunted junction. JoSIM solves the same
RCSJ/capacitor/signed-inductor circuit
in time. Measured source and terminal voltages define the incident, reflected
and transmitted waves. Pump-only traces are subtracted before fitting the
weak probe and idler. Conversion is photon normalized. No coefficients were
fitted to JoSIM response curves.

The port pump is 5 GHz at 1e-10 W incident power; probe power is 1e-15 W.
SQUID/4WM signal and idler are 4 and 6 GHz; rf-SQUID/SNAIL use 2 and 3 GHz.
Loop bias is 0.05 flux quanta. Ports are 50 ohms; each junction and capacitor
has an explicit 10,000 ohm parallel loss resistance. JoSIM uses its phase
formulation, corrector 4 and float64 NPY output. Timestep is 5 fs versus
2.5 fs, output sampling 1 ps, with acquisition at 80–160 ns and a late-window
comparison at 120–160 ns.

Python retains four pump harmonics and three sidebands, compared with six
harmonics and four sidebands plus finer quadrature. Signal/idler differences
must stay below 0.25 dB against JoSIM, while spectral, timestep and acquisition
changes must stay below 0.02 dB. Passive complex scattering error must be
below 0.02 and DC junction-phase error below 1e-4 rad.

| Fixture | Max signal error (dB) | Max idler error (dB) | Idler timestep change (dB) | Max acquisition change (dB) | Max DC phase error (rad) |
| --- | --- | --- | --- | --- | --- |
| SQUID JPA | 0.000029 | 0.000646 | 0.000031 | 0.000683 | 1.27e-08 |
| 4WM JTWPA | 0.000029 | 0.004085 | 0.003818 | 0.000772 | 1.45e-08 |
| rf-SQUID 3WM | 0.000003 | 0.002907 | 0.001996 | 0.002959 | 1.11e-08 |
| SNAIL 3WM | 0.000003 | 0.000405 | 0.000037 | 0.000373 | 2.42e-08 |

A separate SQUID check applies a 0.01-flux-quantum sinusoidal pump through
explicit ideal coupled coils, with zero port-pump power. It passed: the
idler error was 0.000015 dB and DC phase error 7.92e-12 rad.
The prescribed flux coils are an embedding assumption, not an extracted coil mask.

## Rejected attempts and the revision decision

The portable record keeps failed coarse contact comparisons, including the
wider 0.25/0.5 um grids. Intermediate fixed-patch grids were added and the
final adjacent comparisons are listed above; numerical limits were not
relaxed. The earlier 2 um-wire/1 um-contact SNAIL had normalized contact
sensitivity 0.007486, above 0.005. It was not certified. Narrower wiring/contacts
permitted a finer basis within the local size limit, and that new geometry
was independently checked and fully re-extracted. Only its proven identical
capacitor coupons reused electric evidence.

A weaker 4WM test at 1e-11 W pump power had an idler near -124 dB and failed
its timestep comparison because the baseline-subtracted output was too small
for a reliable measurement. That failed record is retained. The final
1e-10 W point resolves the idler and passes both timestep and late-window
checks. Both operating points remain within the AC junction-phase guard.

## KIT calibration-backed check

The two-row KIT example has 270 distributed cells and five period repeats.
Its tapered turns match total biased sheet-squares inductance with relative
error 3.33e-16. Passive energy error is 7.70e-14.
It also runs the existing coupled TLS-loss/amplifier/noise solver; loss tangent
1e-3 and 1% participation are explicitly assumed for this demonstration.
This is a periodic reduction with supplied research C/stub calibration, not
a full serpentine EM extraction or a measurement of dielectric participation.

## Reproduce and interpret

Use the [commented Josephson preparation](../examples/em/layout_network.py),
[KIT example](../examples/em/kit_backannotation.py), and
[native driver](../tools/validate_layout_network.py). The guide gives the exact
grid, geometry, source-impedance and duration choices.
[capture_layout_network_reference.py](../tools/capture_layout_network_reference.py)
checks completed native evidence before writing the portable report. CI
reconstructs the extracted circuits and compares their response with the
stored native measurements without requiring field solvers on every runner.

This is a hybrid component-electric/London-MQS model. Full electric cross-
capacitance, packages, nonlinear film inductance and bias-coil masks are
omitted. The current mesher is bounded to 40,000 filaments and the nodal model
to 200 active terminals; long chains may need a more scalable extraction
method. DC local stability and a small AC phase do not certify dynamical pump
stability. The extracted Josephson response is classical; it does not
automatically inherit the reduced uniform-chain quantum-noise model. KIT
retains its documented forward-wave TLS/noise approximation.

See the [physical background and primary references](references.md),
[numerical-method guide](numerical-methods.md), and the separate
[dedicated JPA benchmark](magnetic-jpa-validation.md).
