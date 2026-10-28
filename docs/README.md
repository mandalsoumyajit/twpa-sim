# TWPA Sim documentation

Use these pages as an engineering workflow: **simulate → lay out → check →
extract → simulate/compare → revise**. The [design-cycle guide](design-workflow.md)
is the organizing route, with a worked JPA iteration and a table showing
which stages are implemented for each family. Start with
[Tools and setup](tools-and-setup.md) to choose the required software for your
path; the Python core, KLayout, JoSIM, openEMS/CSXCAD and superconducting
FastHenry each have a different role.

For the physics behind the stages, read
[What is being simulated?](physical-overview.md). For small electrical
calculations, follow the [HOWTO](../HOWTO.md), then the relevant family guide.
The [examples index](../examples/README.md) maps scripts to the stage they
exercise. The family and validation tables below are references within this
workflow rather than a required front-to-back reading order.

## Work through a design

| Stage | Read / run | What comes next |
| --- | --- | --- |
| Set up the tools | [Required/optional software and runtime checks](tools-and-setup.md) | Install only what the chosen path needs |
| Simulate an electrical candidate | [Design targets and first simulation](design-workflow.md#1-simulate-an-initial-electrical-candidate), [HOWTO](../HOWTO.md) | Revise circuit/drive until the response is useful |
| Generate layout | [Process/layout workflow](layout.md), [worked iteration](design-workflow.md#2-generate-iteration-001) | Keep GDS and its design/process manifest |
| Check geometry/connectivity/DRC | [Check and repair](design-workflow.md#3-check-the-exported-layout-and-fix-geometric-failures) | Regenerate and recheck failures before extraction |
| Extract physical parameters | [Electric component extraction](em-extraction.md), [London network](magnetic-extraction.md) | Refine numerical controls and retain accepted results |
| Back-annotate and compare | [Supported JPA field/circuit run](design-workflow.md#4-extract-back-annotate-and-compare), [JoSIM backend](josim-backend.md) | Check agreement separately from design targets |
| Modify and repeat, or accept | [What to change](design-workflow.md#5-decide-what-to-change), [record the iteration](design-workflow.md#6-keep-a-reproducible-accepted-iteration) | Preserve provenance and repeat affected stages |

Use the [symbol glossary](symbols.md) as a lookup aid, then read the relevant
device guide below. Each guide defines its own equation symbols and states
the physical approximation being solved. The
[numerical-method introduction](numerical-methods.md) explains harmonic
balance, continuation, spectral refinement, stability and noise diagnostics.
[Conventions](conventions.md) is a reference for exact API units, complex
phases and mode ordering, rather than the starting explanation of the physics.

## Device guides

Use this table during the electrical-model stage. Layout generation does not
automatically supply extracted electrical parameters for every family; see
the [stage coverage](design-workflow.md#coverage-by-device-family).

| You want to model | Read | Start with |
| --- | --- | --- |
| Kinetic-inductance line, loading, gradients or compression | [KIT TWPA](kit-twpa.md) | [KIT example](../examples/howto/01_kit.py) |
| Oxide fields, TLS participation and amplifier noise | [TLS pipeline](kit-twpa.md#tls-loss-and-noise) | [Coupled loss example](../examples/howto/02_tls_loss.py) |
| Kerr JPA or effective flux-pumped resonator | [Resonant JPA](resonant-jpa.md) | [Kerr](../examples/howto/03_kerr_jpa.py), [flux](../examples/howto/04_flux_jpa.py) |
| Uniform/RPM junction-chain 4WM | [4WM JTWPA](four-wave-jtwpa.md) | [4WM example](../examples/howto/05_four_wave_jtwpa.py) |
| rf-SQUID/SNAIL 3WM | [Josephson 3WM](three-wave-josephson.md) | [rf-SQUID](../examples/howto/06_rf_squid.py), [SNAIL](../examples/howto/07_snail.py) |
| Quantum scattering, internal baths and added-noise conventions | [Quantum-noise derivation](jtwpa-quantum-noise.md) | Family-specific quantum APIs |
| Physical RCSJ circuits and native/WSL simulation | [JoSIM backend](josim-backend.md) | [Dry-run circuit example](../examples/howto/08_josim.py) |
| Process profiles, mask geometry and KLayout checking | [Layout workflow](layout.md), [process contract](process-contract.md) | [All-family GDS examples](../examples/layout/generate_all.py) |
| JPA capacitor EM extraction and circuit back-annotation | [EM workflow](em-extraction.md) | [openEMS/JoSIM tool](../tools/validate_em_jpa.py) |
| JPA interconnect, kinetic energy and signed inductive coupling | [London-network workflow](magnetic-extraction.md) | [Commented example](../examples/em/london_jpa.py), [full validation tool](../tools/validate_magnetic_jpa.py) |

## Validation and provenance

| Evidence | Report |
| --- | --- |
| Resonant JPA passive/pumped classical response | [JPA/JoSIM](josim-validation.md) |
| 4WM passive scattering, harmonics and depletion | [4WM/JoSIM](jtwpa-josim-validation.md) |
| 4WM Gaussian noise, thermal equilibrium and quantum limits | [4WM noise](jtwpa-quantum-validation.md) |
| rf-SQUID/SNAIL classical response and Gaussian noise | [3WM validation](three-wave-validation.md) |
| KIT source ports, normalization and regression history | [Source audit](source-audit.md) |
| Layout geometry, circuit connectivity and process DRC | [Layout validation](layout-validation.md), [process source review](process-source-review.md) |
| Component-level JPA EM extraction and annotated amplifier response | [EM/JPA validation](em-jpa-validation.md) |
| London-network field convergence, kinetic energy and JoSIM comparison | [Magnetic JPA validation](magnetic-jpa-validation.md) |

The reports state what was tested and its numerical sensitivity. A classical
circuit comparison, quantum consistency test and fabricated-device measurement
are distinct kinds of evidence. The current benchmarks use synthetic or
documented example parameters; they are not a universal process calibration.

Read a report as evidence for a particular circuit and set of numerical
controls. “Residual” measures how well the chosen equations were solved;
“refinement” measures how an observable changes when more numerical detail is
retained. Neither alone establishes physical stability or validates omitted
geometry. A JoSIM comparison checks classical response, while a quantum-limit
check tests the fluctuation model. The [methods guide](numerical-methods.md)
explains why these tests answer different questions.

## Working with the project

- [Units, model boundaries and troubleshooting](conventions.md)
- [Annotated references](references.md) and [BibTeX](../references.bib)
- [Software citation](../CITATION.cff)
- [Contributing and numerical changes](../CONTRIBUTING.md)
- [Josephson implementation scope and roadmap](josephson-roadmap.md)
- [Logo and design provenance](assets/README.md)

The repository uses Markdown documentation so it remains readable directly on
GitHub and in a source checkout. Run `python tools/check_documentation.py`
to check local links/anchors, or add `--run-examples` to execute the HOWTO scripts.
