# Fabrication process JSON contract (v1)

The authoritative semantic validator is
[Process](../src/twpa/layout/process.py). It runs without the optional KLayout
dependency. Load a profile with `load_process('name')` or an explicit JSON path.
Duplicate JSON keys, nonfinite numbers and unknown contract fields are errors.

| Required field | Meaning |
| --- | --- |
| `schema_version` | Integer 1 |
| `id`, `revision`, `description` | Nonempty strings identifying an immutable revision |
| `status` | `research_reference`, `illustrative`, or `local_unqualified`; no signoff status is supported |
| `units` | `um` for all layout coordinates and rule distances |
| `dbu_um` | GDS database unit, positive and finite |
| `grid_um` | Manufacturing-coordinate grid, an integer multiple of DBU |
| `layers` | Role-name mapping to explicit stream pairs and geometry semantics |
| `rules` | Executable rule records with provenance and checking assumptions |
| `sources` | Source ID to human-readable description/locator |
| `parameters` | Numerical geometry/material reference records, with explicit units and sources |
| `limitations` | Nonempty list of unresolved calibration, stack and checking boundaries |

Every layer record contains `gds: [layer, datatype]`, `purpose`, `semantics`
and `source`. Stream mappings must be unique integer pairs. A missing datatype
cannot be guessed as zero. `semantics` is one of `conductor`, `opening`,
`junction`, `capacitor`, `inductor`, `boundary`. `DIE` is required. Generator
roles are documented in the shipped profiles. This contract describes retained
conductors, openings and annotations; photoresist mask tone/inversion is a
separate process-adapter responsibility.

Every rule contains `id`, `kind`, ordered `layers`, `value_um`, `source`, and
`status` (`reference_geometry` or `engineering_assumption`). Rule IDs must be
unique and all layer/source references must resolve. Unsupported semantics
must stay in `limitations`, not be encoded as a rule that silently passes.

| `kind` | Layers, in order | Exact check |
| --- | --- | --- |
| `width` | One layer | KLayout Euclidean minimum width |
| `space` | One layer | KLayout minimum spacing, including notches within a merged polygon |
| `grid` | One layer | Every unmerged polygon vertex on the specified coordinate grid |
| `inside` | Inner, outer | Boolean inner-minus-outer is empty; `value_um=0` |
| `enclosure` | Inner, outer | Containment **and** the positive minimum enclosure margin |

Distances are integer multiples of DBU. Width/space may specify optional
`ignore_angle_deg` in (0,90]; edges meeting at angles at or above that value
are excluded by KLayout's angle filter. Default 90 degrees. The KIT profile
explicitly uses 80 degrees because grid-snapped rotated rectangles acquire
apparent slightly acute corners. Parallel opposing edges and real narrow
width/spacing defects are still checked. This is an explicit research
checking convention, not a vendor exemption.

Parameter records contain `value`, `unit`, `source`, `status`. Examples are
`sheet_inductance` in `pH/sq`, `critical_current_density` in `uA/um2`, and
specific capacitance in `fF/um2`. `Process.parameter(name, unit)` rejects
unit mismatches. The code applies no implicit process bias, rounding,
statistical corner or density recalibration. The constrained electrical
sizing helper explicitly rounds dimensions to the declared grid and checks
the resulting electrical error. User-supplied geometric coordinates must
already lie on-grid. Curved KIT vertices are explicitly snapped to DBU;
the current curved generator requires `grid_um == dbu_um`.

## Design and export contracts

Design JSON contains exactly `family` and `geometry`. Families and supported
dataclass fields are in [generate.py](../src/twpa/layout/generate.py).
Unknown fields and inconsistent dimensions fail before export. Geometric
constraints define the support envelope of these initial primitives; changing
a process's minimum rules does not automatically resynthesize legal devices.

Every `.manifest.json` carries the complete process snapshot, typed design,
canonical SHA-256 hashes, exported GDS hash, terminal coordinates, component
intent, area-based estimates and extraction limitations. Coordinate and
marker units are um/um²; `Ic_A`, `Cj_F` and `C_F` are SI. Process/design hashes
must match the stored objects. The GDS hash is informative so manual edits
can be checked against unchanged intent; `gds_matches_manifest` records it.

Validation uses serialized polygons. It compares actual JJ/CAP/IND geometry
with component intent, identifies conductor islands and merges M1/M2 only
through VIA openings. A JJ joins an M1 terminal and an M2 terminal as a circuit
edge, not a galvanic short. CAP joins its signal M2 electrode to ground M1 as
a circuit edge. IND cuts M1 in the extraction view only, exposing the two
inductor terminals. Expected net assignments from pins/components must agree
with the extracted conductor partition. Extra floating islands, intended-net
opens and shorts between distinct intended nets fail.

Reports provide rule-by-rule counts, geometry violations, physical/intended
net counts, process qualification, hashes and KLayout version. `.lyrdb`
contains spatial markers for geometric rules and text items for topology
errors. `drc_executed` distinguishes full checking from geometry-only mode.

## Extending a real process

Create a separately versioned local profile from reviewed source documents.
Resolve stream datatypes, polygon tone, legal contact truth tables, junction
fabricated-area bias, grids, enclosure contexts and units. Add new generators
or adapters when the layer roles or stack differ; renaming a layer map cannot
turn this overlap example into an anodized-junction process. Add executable
conditional checks only after their semantics are resolved. Keep the
unimplemented/manual-rule inventory explicit and qualify any claim of complete
coverage independently of the current research profiles.
