# TWPA Sim documentation

Start with the [HOWTO](../HOWTO.md), which links eight commented, runnable
examples. Use [conventions](conventions.md) to compare APIs safely, then
follow a device guide for equations and numerical controls.

## Device guides

| You want to model | Read | Start with |
| --- | --- | --- |
| Kinetic-inductance line, loading, gradients or compression | [KIT TWPA](kit-twpa.md) | [KIT example](../examples/howto/01_kit.py) |
| Oxide fields, TLS participation and amplifier noise | [TLS pipeline](kit-twpa.md#tls-loss-and-noise) | [Coupled loss example](../examples/howto/02_tls_loss.py) |
| Kerr JPA or effective flux-pumped resonator | [Resonant JPA](resonant-jpa.md) | [Kerr](../examples/howto/03_kerr_jpa.py), [flux](../examples/howto/04_flux_jpa.py) |
| Uniform/RPM junction-chain 4WM | [4WM JTWPA](four-wave-jtwpa.md) | [4WM example](../examples/howto/05_four_wave_jtwpa.py) |
| rf-SQUID/SNAIL 3WM | [Josephson 3WM](three-wave-josephson.md) | [rf-SQUID](../examples/howto/06_rf_squid.py), [SNAIL](../examples/howto/07_snail.py) |
| Quantum scattering, internal baths and added-noise conventions | [Quantum-noise derivation](jtwpa-quantum-noise.md) | Family-specific quantum APIs |
| Physical RCSJ circuits and native/WSL simulation | [JoSIM backend](josim-backend.md) | [Dry-run circuit example](../examples/howto/08_josim.py) |

## Validation and provenance

| Evidence | Report |
| --- | --- |
| Resonant JPA passive/pumped classical response | [JPA/JoSIM](josim-validation.md) |
| 4WM passive scattering, harmonics and depletion | [4WM/JoSIM](jtwpa-josim-validation.md) |
| 4WM Gaussian noise, thermal equilibrium and quantum limits | [4WM noise](jtwpa-quantum-validation.md) |
| rf-SQUID/SNAIL classical response and Gaussian noise | [3WM validation](three-wave-validation.md) |
| KIT source ports, normalization and regression history | [Source audit](source-audit.md) |

The reports state what was tested and its numerical sensitivity. A classical
circuit comparison, quantum consistency test and fabricated-device measurement
are distinct kinds of evidence. The current benchmarks use synthetic or
documented example parameters; they are not a universal process calibration.

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
