# JPA EM extraction and circuit back-annotation

`twpa.em` starts the layout → EM → circuit workflow with the single-junction
resonant JPA. It extracts a **MIM capacitor coupon from the exported GDS**,
solves it with openEMS, checks numerical refinements, and replaces the nominal
shunt capacitance in a physical RCSJ circuit. The resulting circuit feeds the
existing reduced JPA model and independent JoSIM transient calculations.

This page describes the initial component-level electric extraction. The
[London-network extension](magnetic-extraction.md) adds finite interconnect,
ground return, kinetic energy and signed inductive coupling. Package and
bias geometry remain outside those benchmarks. Junction `Ic` and `Cj` come from actual
GDS junction-window area and the explicit process densities. PEC extraction
does not determine superconducting loss or nonlinear junction physics.

## Physical stack

The [illustrative EM stack](../src/twpa/em/stacks/mim_illustrative.json) defines
the dielectric thickness/permittivity, substrate thickness/permittivity and
PEC-sheet electrode model. It is a separate physical contract from the
[mask process profile](process-contract.md). The default dielectric corresponds
approximately to 1 fF/um² in the infinite parallel-plate limit. This is an
illustrative stack, not a calibrated foundry material model.

`EMStack` rejects unsupported electrode models, qualification claims and
nonpositive/nonfinite physical inputs. Supply a local stack JSON to change
materials deliberately. There is no implicit inference of dielectric
thickness from nominal circuit capacitance.

## Geometry and reference planes

`capacitor_job(gds, stack)` first requires the original GDS/manifest hash match
and a clean geometry/connectivity/DRC check. It selects the JPA's single CAP
component, crops actual M1/M2 polygons around that component, and exports
hole-free polygon pieces to CSXCAD. The default crop margin is 1 um beyond
the CAP annotation. Ports and electrode coordinates are in the capacitor
coupon, not at the full amplifier's external RF terminal.

In the default layout, the annotation is 12 × 12 um², while actual metal
overlap includes electrode enclosure and a clipped feed neck. The extractor
therefore does not simply substitute `annotation_area * specific_capacitance`.
The coupon omits the rest of the leads; its crop and footprint are included
in the extraction artifact for review. It requires one connected electrode
per metal and rejects unsupported families and ambiguous crops.

openEMS excites a small differential lumped port between the two PEC sheets.
All polygon vertices, port edges and dielectric interfaces are placed on
the mesh. A uniform lateral grid and explicit dielectric subdivisions define
the baseline, with nonuniform vertical air/substrate planes. Boundaries are
MUR. The default port resistance is a **numerical extraction termination**;
it is not the amplifier's external impedance or internal loss.

The port admittance uses openEMS's `exp(+i omega t)` convention. The extracted
real capacitance is obtained from `Im(Y)/(2*pi*f)` and a parallel-C fit across
the frequency band. `twpa` circuit calculations use `exp(-i omega t)`; a real
capacitance transfers without a sign ambiguity. Small residual real
admittance is recorded as a numerical/domain diagnostic and is not converted
into a material-loss resistor.

## Numerical acceptance and provenance

A single successful process exit is insufficient. Each extraction requires:

- Positive finite capacitance and a maximum complex-admittance fit residual
  below 2% across the supplied band.
- Voltage/current port tails below 0.1% of their peaks.
- Negative real admittance below 1% of total admittance magnitude.
- Separate comparisons of lateral/vertical mesh, domain margin, and duration.

`certify(primary, comparisons)` checks physical geometry, source, stack and
solver identity, and requires all three comparison groups. The primary must
use the finest mesh, largest domain and longest duration in the study.
Default acceptance permits at most 2% change in capacitance at every sampled
frequency; it is a convergence sensitivity gate, not a universal error bound.
The reference run reaches its fixed timestep limit with residual global
static-field energy. Acceptance relies on converged **port** tails and the
independent duration comparison; it does not claim that the global -50 dB
FDTD energy criterion was met.

Each fresh run directory retains `job.json`, the exact worker source,
`model.xml`, voltage/current probes, solver logs, runtime metadata and
`extraction.json`. Jobs, worker source, layout, process, design and stack
are hashed. Worker failures/timeouts preserve diagnostics; existing run
directories cannot be overwritten. Extraction artifacts have a canonical
integrity hash. Stale GDS or manifest hashes reject back-annotation.

## Run the complete workflow

Install `twpa-sim[layout]` in the main Python environment. Install openEMS and
CSXCAD using their [official instructions](https://docs.openems.de/); their
compiled Python modules may need a separate Python version. The worker runs
in that interpreter, so it need not contain `twpa`, KLayout or SciPy.

The [validation tool](../tools/validate_em_jpa.py) generates the JPA, runs four
EM cases, certifies extraction, writes the annotated circuit, and compares
passive reflection, signal gain and idler conversion against JoSIM:

```powershell
python tools/validate_em_jpa.py --output validation/runs/my-em-jpa --openems-python C:/path/to/openems/python.exe --openems-dll-directory C:/path/to/openEMS --executable C:/path/to/josim-cli.exe
```

On Linux, omit `--openems-dll-directory`. To run JoSIM in WSL, add `--wsl Ubuntu`
and use a Linux executable path. Fork-specific `--corrector` is optional and
checked against the selected executable. `--revision` records a build source
identifier; record the exact binary hash for archival validation as well.

Use `--extract-only` to stop after EM/back-annotation. To reuse a certified
result, pass `--gds path/to/jpa.gds --extraction path/to/extraction-certified.json`
and a fresh `--output` directory. The full validation tool saves compact
measurements and all circuit inputs/logs; completed JoSIM CSV waveforms are
removed after loading unless `--retain-waveforms` is supplied. Generated
artifacts belong under ignored `validation/runs/`.

## Python interface

```python
from twpa.em import load_stack, capacitor_job, OpenEMS

job = capacitor_job('jpa.gds', load_stack(), mesh_um=0.5, z_cells=3)
solver = OpenEMS('/path/to/openems/python')
result = solver.run(job, 'validation/runs/capacitor-primary')
# Run separate mesh/domain/duration comparisons before certify(...).
```

After certification:

```python
from twpa.em import back_annotate_jpa, load_circuit

circuit, record = back_annotate_jpa(
    'jpa.gds', 'extraction-certified.json',
    port_impedance=5000, internal_resistance=1e7,
    output='backannotation.json')
assert load_circuit('backannotation.json') == circuit
resonator = circuit.reduced_model()
```

The annotated `ParallelJPA` works with `jpa` response functions and
`jpa_netlist`/`run_circuit`. Junction capacitance is added exactly once,
separately from the EM shunt. External impedance, internal resistance and
junction resistance remain explicit user assumptions with per-parameter
origins recorded in `backannotation.json`. The reduced resonance must lie
inside the extraction frequency band.

See the [initial capacitor validation](em-jpa-validation.md) and the newer
[London-network validation](magnetic-jpa-validation.md), which retains finite
interconnect, return, contacts, kinetic energy and inductive coupling. Electric
stray capacitance beyond the MIM coupon, packages, dielectric/material loss,
bias structures and SQUID/SNAIL/JTWPA back-annotation remain separate work.
Neither benchmark establishes a fabricated amplifier's measured performance.
