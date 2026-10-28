# KIT implementation provenance

The KIT routines were ported from earlier MATLAB and Python design scripts.
This audit maps those sources to the package and records numerical and
normalization changes. It preserves implementation history, rather than
serving as an independent experimental or literature validation.

| Original | Extracted implementation |
| --- | --- |
| `sim_tls_and_spine.py` | ABCD, complex Bloch propagation, calibrated Cu/Cl, 925-period circuit, constant-pump lossy gain and Langevin noise |
| `generate_reorganized_figures.py` | `compute_2d_laplace_fields`: conformal oxide geometry, finite-volume stencil, boundaries and field energy |
| `sim_alvin_twpa.py` | Complete nonlinear three-current RHS with pump depletion and Kerr phases |
| `sim_detailed_physics.py` | Same nonlinear RHS, compression sweep and power conventions |
| `sim_gradients_exact.py` | WKB phase-integral gradient model and spatial profiles |
| `generate_presentation_figures.py` | Alternative symmetric cell topology |
| `kitwpa_nycreates_v1.m` | Four historical electrical presets |
| `KI_TWPA_Design_Report.md` | Final device counts and pipeline prescription |
| `generate_v2_gds.py` | Parameterized straight, U-turn and launch primitives in `twpa.layout.kit`; [geometry/provenance audit](layout-validation.md) |

The coupled TLS solver is present in the earlier sources. It uses a constant
pump and lossy signal/idler Bloch parameters in a rotating-frame transfer and
Green-function noise integral. `constant_langevin` preserves this; `coupled`
adds optional pump attenuation and thermal baths.

The original figure script solves a Laplace field image, but plots tabulated
participation/capacitance curves. `oxide_design` implements the report's energy
ratio prescription to connect actual field extraction to the amplifier. This
is new wiring of ported components, not a claim that every original plotted
curve was automatically computed through all three stages.

Deliberate differences:

* Explicit Kerr factor preserves the source's factor-of-two discrepancy.
* Source Cu/Cl are total capacitance; per-finger Cell values divide them by two.
* Continuous Bloch branching replaces the fixed 8.46 GHz threshold; the full
  trace replaces the A-only shortcut used in older dispersion scripts.
* Scaled matrix powers prevent overflow. Zero-growth gain avoids division by
  zero. Constant bath covariance uses exact matrix exponentials rather than
  finite trapezoidal integration.
* Gradient phase is integrated directly. The source spiral expression labeled
  an exact antiderivative is not the integral of its supplied coordinate.
* Compression reports NaN without a bracketed crossing rather than silently
  returning an endpoint through extrapolating/clamping interpolation.
* Import-time plotting, GDS work and personal absolute output paths are removed.

`tools/capture_reference.py` executes only original numerical definitions and
the TLS sweep, skipping plot/GDS side effects. Its fixture records source
SHA256 hashes, all six TLS cases, lossy/lossless Bloch samples, original 2.5 nm
Laplace energy/participation/potential probes, and nonlinear current outputs.
Tests check source gain at rtol 1e-11 and bath noise at 2e-6 (the old quadrature
is approximate), plus Bloch, fields, nonlinear currents and coupled ODE results.

Additional tests cover analytic attenuation and noise limits, depletion and
compression, zero-gradient equivalence, oxide-to-amplifier execution, custom
coefficients, topology, scattering and convergence. Source report performance
claims and hardcoded oxide tables are not independent numerical verification.
