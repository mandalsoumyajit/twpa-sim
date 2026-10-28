# Examples: what to run and what it means

Run scripts from the repository root after installing `twpa-sim`. Begin with
the [physical introduction](../docs/physical-overview.md), then follow the
[HOWTO](../HOWTO.md). The numbered scripts below are the smallest learning
cases and are exercised by the documentation check. Their opening text and
comments explain the device, parameter choices and output quantities.

All electrical values are illustrative unless a guide explicitly identifies
a literature preset. A short example with modest gain is not a failed attempt
to reproduce a published high-gain amplifier. It is intended to expose the
calculation and make parameter changes inexpensive.

## Guided electrical examples

| Script | Circuit and calculation | First outputs to inspect |
| --- | --- | --- |
| [01_kit.py](howto/01_kit.py) | Biased kinetic-inductance line; undepleted versus finite-signal propagation | Weak/finite gain, remaining pump current, passive power balance |
| [02_tls_loss.py](howto/02_tls_loss.py) | Surface-oxide electric participation feeding a coupled loss/amplifier model | Participation, capacitance ratio, gain and input-referred added noise |
| [03_kerr_jpa.py](howto/03_kerr_jpa.py) | Effective nonlinear resonator, pumped through its measurement port | Stable pump occupation, reflected signal gain and added noise |
| [04_flux_jpa.py](howto/04_flux_jpa.py) | Effective flux modulation at twice resonance; coincident signal and idler | Two quadrature gains and their noise relative to vacuum 0.5 |
| [05_four_wave_jtwpa.py](howto/05_four_wave_jtwpa.py) | Finite unbiased junction chain with resistive losses | Transmitted gain, idler photon conversion, noise bound and quantum consistency |
| [06_rf_squid.py](howto/06_rf_squid.py) | Magnetically biased one-junction loops in a line | Potential derivatives, average-phase shift, transmitted gain and noise |
| [07_snail.py](howto/07_snail.py) | Asymmetric small/large-junction loops in a line | Equilibrium phase, potential derivatives, gain and noise |
| [08_josim.py](howto/08_josim.py) | Physical junction netlist; optional classical transient run | Saved components/connections, then voltage and phase traces if executed |

Gain in dB is `10*log10(G)`, where `G` is output signal power divided by
incident signal power. Josephson idler conversion is a photon-flux ratio;
photon flux is power divided by photon energy. Added noise is referred to
the input and excludes the original signal input's noise. The
[HOWTO output guide](../HOWTO.md#reading-the-output) explains these definitions
and the finite-gain quantum limit.

## Geometry and field extraction

Install `twpa-sim[layout]` for these examples:

```powershell
python examples/layout/generate_all.py
python examples/em/london_jpa.py --output validation/runs/my-jpa
```

[generate_all.py](layout/generate_all.py) reads process/design JSON files,
writes GDS geometry and manifests, and checks connectivity and design rules.
It does not simulate amplifier gain. A full KLayout executable can additionally
run the independent native design-rule deck; see [the layout guide](../docs/layout.md).

[london_jpa.py](em/london_jpa.py) prepares the finite metal-current network of
a resonant JPA with a physical pickup loop. `--bare` instead uses its original
direct electrical port. It prints conductor area, current-filament count and
the oriented terminal pairs, and writes a FastHenry input deck. It does not
launch a field solve. Continue with the [extraction guide](../docs/magnetic-extraction.md)
for material inputs, numerical refinement, electric capacitance extraction
and circuit back-annotation. A locally generated deck is not yet an accepted
inductance matrix.

## Longer exploratory scripts

These scripts explore more cases or larger circuits. Read the related guide
and complete the numbered example first; runs can be more expensive.

| Script | Purpose | Background |
| --- | --- | --- |
| [compare.py](compare.py) | Compare historical KIT electrical presets over signal frequency | [KIT model](../docs/kit-twpa.md) and [source audit](../docs/source-audit.md) |
| [depletion_and_gradients.py](depletion_and_gradients.py) | Explore finite-signal propagation with assumed spatial inductance profiles | [KIT model](../docs/kit-twpa.md) |
| [electrostatic_input.py](electrostatic_input.py) | Demonstrate a small dielectric cross-section boundary-value calculation | [KIT/TLS fields](../docs/kit-twpa.md#tls-loss-and-noise) |
| [oxide_pipeline.py](oxide_pipeline.py) | Feed an oxide electric-energy calculation into a larger lossy KIT example | [TLS model](../docs/kit-twpa.md#tls-loss-and-noise) |
| [tls_amplifier.py](tls_amplifier.py) | Study the coupled loss/amplifier calculation with prescribed participation | [TLS model](../docs/kit-twpa.md#tls-loss-and-noise) |
| [resonant_jpa.py](resonant_jpa.py) | Sweep Kerr gain and effective flux quadratures, with optional plotting | [JPA model](../docs/resonant-jpa.md) |
| [four_wave_jtwpa.py](four_wave_jtwpa.py) | Compare uniform and resonantly loaded finite junction chains | [4WM model](../docs/four-wave-jtwpa.md) |
| [quantum_jtwpa.py](quantum_jtwpa.py) | Compare vacuum and thermal internal baths in a larger lossy 4WM chain | [Noise derivation](../docs/jtwpa-quantum-noise.md) |
| [three_wave_josephson.py](three_wave_josephson.py) | Compare synthetic rf-SQUID and SNAIL line responses | [3WM model](../docs/three-wave-josephson.md) |

The [annotated bibliography](../docs/references.md) distinguishes physical
background papers from numerical implementation and validation evidence.
Use the [methods guide](../docs/numerical-methods.md) to refine a run before
interpreting a parameter trend quantitatively.
