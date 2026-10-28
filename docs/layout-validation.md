# Layout validation

Validation performed on 2026-10-27 using standalone KLayout **0.30.12** and
the installed native KLayout **0.30.9** Ruby engine. The profiles are a research
KIT reference and an illustrative JJ stack; passing these checks is not
foundry signoff or an EM validation of amplifier performance.

## Original KIT primitive regression

[capture_layout_reference.py](../tools/capture_layout_reference.py) executes
only the original `get_v2_primitives` definition and compares its regions
with the parameterized default construction. Straight supercell, right and
left U-turns, and both launches all have **zero XOR area** at 1 nm DBU.
The [record](../validation/layout-reference-primitives.json) pins the source
SHA-256, areas, three islands per primitive and solved turn maximum width
1.1571445906404483 um. Regression tests preserve the reference areas/counts.
This establishes primitive fidelity, not identity of the entire exported mask:
landing-pad placement is generalized and the other original routing variants
are not yet ported.

## Full reference size

The default 16-row, 55-period-per-row KIT has **49,950 distributed cells**.
Its exported GDS passes Python geometry/connectivity and every profile DRC
rule. Native Ruby DRC also reports **zero violations**. It has the intended
three electrical conductor islands, with no unassigned islands. The two
ground combs remain distinct and need package connections.

On this workstation, generation plus Python reimport/geometry/DRC took about
159 seconds; native DRC was run separately. Small fixtures are used in CI
to keep regressions practical. Full GDS and spatial marker files are local
regenerable artifacts under ignored `validation/runs/layout`.

## Device fixtures and planted failures

[test_layout.py](../tests/test_layout.py) checks all six device families,
changed KIT dimensions/loading/odd row count, changed SNAIL junction sizes and
array count, and electrical sizing consistency. Traveling-wave capacitor
areas sum to `cells * ground_capacitance`, using endpoint half-shunts.
Every positive fixture has consistent intended and physical net partitions.

The checker must reject deliberately planted:

- Signal-to-ground shorts and breaks in intended nets.
- Missing junction markers, contacts or inductor annotations.
- Broken conductor inside an inductor extraction window.
- Too-narrow metal and too-close spacing.
- Contacts outside their enclosing metals.
- Off-grid vertices, geometry outside the die and unknown stream layers.
- Altered process/design hashes and malformed process contracts.
- Incompatible junction current/capacitance targets during electrical sizing.

Native Ruby DRC separately checks the positive family fixtures and planted
width, space, enclosure and grid failures. The Python validator's topology
checks are necessary in addition to native DRC: legal widths alone cannot
establish an intended circuit. A zero native process exit code is never used
as the pass criterion; the wrapper reads `.lyrdb` items.

The compact [all-family evidence record](../validation/layout-validation.json)
contains process/design/GDS hashes, rule counts, qualifications, component
counts and extracted net counts. Regenerate it with:

```powershell
python examples/layout/generate_all.py --klayout klayout --summary validation/layout-validation.json
```

Use the full executable path on Windows. Add `--full-kit` to replace the small
KIT fixture with the full reference. A [separate full-reference record](../validation/layout-full-reference.json)
summarizes the full-size run. Run unit tests with:

```powershell
python -m unittest discover -s tests -p test_layout.py -v
```

Set `TWPA_KLAYOUT_EXECUTABLE` to enable the native integration test. CI has a
dedicated layout job that installs both the standalone module and a full
KLayout package, runs those tests, then generates/checks the JSON examples.
The core simulation tests can run without layout dependencies.

## Remaining implementation

The initial workflow deliberately exposes these boundaries:

- KIT folded/four-pass routing and configurable package-wide ground stitching.
- Foundry-specific JJ/barrier/contact/capacitor primitives and conditional DRC.
- Mask-tone conversion, junction etch bias, and process-qualified corners.
- Calibrated SQUID/SNAIL loop and TWPA cell extraction/back-annotation.
- RPM resonators, electric coupling networks, bias coils/tees and on-chip loss elements.
- Package extraction and electric stray capacitance outside the MIM coupon.

The [JPA London-network implementation](magnetic-extraction.md) now covers
the capacitor coupon plus finite interconnect, return, contact, kinetic energy
and signed inductive coupling extraction. The new pickup-loop geometry is
included in the all-family geometry and native DRC example campaign.

The earlier JoSIM reports apply to their numerical circuits. The separate
[London-network report](magnetic-jpa-validation.md) checks the extracted
full-sine circuit against JoSIM. Omitted electric/package effects and other
device families remain outside that benchmark.
