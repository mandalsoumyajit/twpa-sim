# Josephson model implementation scope

Status: resonant JPA models and finite-circuit 4WM JTWPA models implemented
and compared with JoSIM. Uniform/RPM chains include multiharmonic mixing and
classical pump depletion. See the [4WM guide](four-wave-jtwpa.md) and
[circuit validation](jtwpa-josim-validation.md). The targets below preserve
the broader scope. Gaussian 4WM quantum noise with internal resistor baths
is now implemented; see [quantum noise](jtwpa-quantum-noise.md).
Flux-biased rf-SQUID/SNAIL 3WM potentials, DC rectification, classical
depletion and pump-only Gaussian quantum noise are implemented; see the
[3WM guide](three-wave-josephson.md) and [validation](three-wave-validation.md).
Updated 2026-10-27. See [implemented JPA models](resonant-jpa.md).

## Targets and order

Implement a conventional resonant JPA first, then a junction-chain 4WM JTWPA,
then flux-biased SQUID/SNAIL 3WM traveling-wave devices. JoSIM is the optional
classical circuit-simulation backend for each Josephson design. Keep the
existing 19 KIT/circuit regression tests as compatibility gates.

“Plain JPA” here means a single resonant Josephson parametric amplifier in
reflection. Include current-pumped Kerr operation and an effective flux-pumped
degenerate parametric oscillator. A Josephson parametric converter, multimode
resonator, or impedance-engineered broadband JPA is a subsequent extension.

## Existing code: reuse and gaps

Reuse scaled matrix cascades, finite-network scattering, loading sequences,
ODE utilities, pump sweeps and bath-covariance machinery. Preserve current
public APIs through adapters rather than changing the meanings of Cell fields.

The existing Cell is an IDC/stub ladder, not an arbitrary JJ circuit. Introduce
a circuit-element protocol returning normalized ABCD matrices versus frequency.
The current nonlinear/gain functions impose 3WM frequency conservation and a
particular current-envelope interaction. A custom Kerr coefficient alone does
not make them 4WM solvers. Extract a common transfer/noise integrator accepting
explicit coupling, detuning and attenuation rather than internally deriving
them from the existing Bias law.

Keep frequency [Hz] distinct from angular frequency/rates [rad/s], cell index
from length, peak current from photon-flux amplitudes, and attenuation from
nonlinear phase shift. Existing circuits use exp(+i omega t); input-output
equations below use exp(-i omega t). Supply explicit convention conversions.

## 1. Resonant JPA

Suggested module: `twpa.jpa`.

Inputs: resonance frequency, Kerr K, external/internal linewidths, bath
temperatures, pump frequency/power/phase, and signal detuning. Use cavity
occupation |a|^2 and incident photon flux |s_in|^2, with
P_in = hbar*omega_p*|s_in|^2. Require linewidth units to be explicit.

For the Kerr model, use Delta = omega_p - omega_0 and

    da/dt = (i Delta - kappa/2) a - i K |a|^2 a
            + sqrt(kappa_ext) s_in

where kappa = kappa_ext + kappa_int. Solve every nonnegative real root of

    n * [(kappa/2)^2 + (Delta - K*n)^2]
        = kappa_ext * P_in/(hbar*omega_p)

Reconstruct the complex pump amplitude, linearize about each solution, and
classify stability from drift eigenvalues. Report stable/unstable branches;
do not silently choose a root inside a bistable region. Supply continuation
when a sweep is intended to follow a branch. Compute the reflected signal,
conjugate-idler conversion, bandwidth and internal-loss bath contributions.

The effective flux-pumped model takes a complex pair-drive lambda and rotates
at half the physical pump frequency. Return quadrature gains and squeezing
at degeneracy, rather than treating coincident signal/idler as independent
vacuum inputs. The stability boundary for that linear oscillator is
|lambda|^2 = (kappa/2)^2 + Delta_eff^2. Physical SQUID-flux-to-lambda extraction
is a separate circuit adapter, not an assumed universal conversion.

Acceptance: zero-pump passive reflection and linewidth; all pump roots satisfy
the residual; stability/threshold checks; lossless commutator preservation;
finite-gain phase-preserving noise limit away from degeneracy; degenerate
quadrature uncertainty; thermal-loss behavior; time-domain compression and
hysteresis consistent with the selected branch. Benchmark against the
[Eichler–Wallraff input-output treatment](https://link.springer.com/article/10.1140/epjqt2).

## 2. Junction-chain JTWPA, 4WM

Suggested modules: `twpa.josephson`, `twpa.mixing` and `twpa.jtwpa`.

Inputs: junction Ic, Cj, subgap damping/model parameters, operating phase/flux,
ground capacitance, stray inductance, cell length, number of cells, ports and
optional capacitively coupled LC phase-matching resonators. Start from the
Josephson potential/current-phase relation and a verified equilibrium branch.
Linearized Lj = (Phi0/2pi)/(Ic*cos(phi_bias)); include Cj and damping rather
than substituting Lj into the IDC finger circuit. Initially support stable,
positive differential-inductance operating points.

Implement 4WM with omega_i = 2*omega_p - omega_s and linear mismatch
2*k_p - k_s - k_i. Derive device-specific mixing and self/cross-phase shifts
in one documented normalization; the pump-mediated pair coupling is quadratic
in pump amplitude. Validate an undepleted solution, then add consistent
pump-depleted propagation, loss and quantum bath covariance. Keep dissipative
attenuation distinct from the phase shifts sometimes also denoted alpha.

Acceptance: unpumped dispersion versus direct circuit analysis and JoSIM;
junction plasma-frequency limit; pump-off passivity; 4WM frequency and
Manley–Rowe balances in the declared amplitude normalization; mismatch and
gain versus the [O’Brien et al. phase-matching model](https://arxiv.org/html/1406.2346v1);
depletion, compression and harmonics versus transient simulation. Reproduce
the stated circuit and a selected gain curve before adding fabrication spread.
Generic gain/noise tests should not substitute for this circuit benchmark.

## 3. SQUID/SNAIL 3WM JTWPA

Represent the flux-dependent potential with junction asymmetry, loop
inductance/screening where relevant, and a chosen stable equilibrium. Extract
quadratic, cubic and quartic derivatives, with a documented conversion into
mode/flux or current amplitudes. Retain residual Kerr effects and include
loading or resonators that control pump harmonics. Validate first in the weak
drive limit; do not equate a numerically small bare quartic coefficient with
a universally Kerr-free effective device.

Acceptance: equilibrium residual and positive curvature; periodicity under a
flux-quantum shift; symmetry-forbidden cubic coefficient where applicable;
harmonic content and small-signal gain versus JoSIM. Use
[Roudsari et al.](https://arxiv.org/abs/2209.07551) as a 3WM benchmark candidate.

## 4. JoSIM backend and measurement extraction

Suggested modules: `twpa.backends.josim`, `twpa.measurement`.

Local checkout inspected: `C:\dev\JoSIM`, commit
`5f33997c7c8a7317e98e3aa8c49ce4503375d9a2`, fork
`mandalsoumyajit/JoSIM`. It contains an enhanced RCSJ transient simulator,
SPICE-style parameterized decks and WSL build scripts. Its WSL CLI build has
now been used for resonant JPA validation. See [backend configuration](josim-backend.md)
and [measured comparisons](josim-validation.md). No developer checkout path
is hardcoded in the backend.

Backend contract:

1. Configure executable and optional WSL launcher explicitly; detect version
   and capabilities. Use argument arrays, run-directory isolation, timeout,
   return-code checking, captured logs and reproducible random seeds.
2. Generate RCSJ junction-chain, resonator, SQUID/SNAIL and matched-port decks
   from the same physical parameters used by reduced models. Avoid counting
   Cj twice if it is already present in the junction model. Flux bias/pumping
   needs explicit loop inductors and mutual inductances, not a fictitious JJ
   flux terminal. Pin boundary conditions and model parameters.
3. Start with upstream-compatible CSV and transient controls. Probe before
   using this local fork's NPY/binary output, checkpointing or corrector
   options. Checkpoint reuse must match the netlist, bias and pump history.
4. Establish pump-only steady state, ramp a weak probe, discard settling and
   fit complex tones coherently. Subtract pump-only leakage where needed.
   Compute incident/reflected port waves from voltage and current with a
   consistent reference plane; output/input voltage ratios alone are not
   transducer gain. Extract the generated idler, harmonics and power balance.
5. Sweep probe power for compression; sweep pump and initial conditions for
   bistability. Demonstrate timestep, settling-time and acquisition-window
   convergence; resolve the highest retained harmonic/plasma oscillation,
   not just the signal frequency. Check spatial/device-count scaling before
   launching full thousand-junction sweeps.
6. Preserve netlist, simulator revision, parameters, source phases, solver
   controls, measurement windows and extracted results per run. Add small
   optional integration tests, skipped when JoSIM is unavailable, and frozen
   fixtures with provenance. Keep JoSIM an external dependency; do not vendor
   its implementation or assume the fork's optional features exist upstream.

[JoSIM supports transient JJ simulation](https://joeydelp.github.io/JoSIM/)
and [classical resistor thermal noise](https://joeydelp.github.io/JoSIM/tech_disc/).
That noise model is not vacuum noise or a prediction of quantum-limited added
noise. Use JoSIM to validate classical dynamics and calibrate reduced models;
retain quantum input-output/Langevin models for squeezing and added photons.

## Deliverable milestones

| Milestone | Reviewable result and exit gate |
| --- | --- |
| M0 | Circuit/mixing interfaces and normalization document; existing 19 tests unchanged |
| M1 | Kerr and effective flux-pumped JPA, roots/stability/scattering/noise; analytic limits pass |
| M2 | JoSIM runner, RCSJ/JPA decks and coherent port-wave extraction; converged JPA comparison |
| M3 | Junction-chain/RPM linear circuits and 4WM reduced solver; published-circuit and JoSIM checks |
| M4 | JTWPA depletion/compression/loss/noise; balance and convergence gates |
| M5 | SQUID/SNAIL potential extraction and 3WM adapter; symmetry and transient benchmarks |

The first JPA release includes stable pump roots, effective flux pumping,
sideband scattering, quantum noise, quadratures and classical probe transients.
JoSIM validation now covers passive reflection, Kerr gain/idler/compression,
and flux-pump sidebands and quadratures. A dedicated resonant compression
search and dynamic hysteresis sweeps remain future work. The JTWPA release
implements M3's finite-circuit uniform/RPM dispersion and gain, plus M4's
classical depletion, resistive loss and Gaussian quantum bath covariance.
It uses `twpa.jtwpa` harmonic balance
rather than the originally suggested continuum transfer refactor. The paper
circuit preset is supplied, but full published gain-curve reproduction,
automated compression searches, arbitrary fabrication patterns and dynamical
stability analysis remain future work. These APIs
are separate and preserve the existing KIT and JPA interfaces.

M5 now includes a unique-branch rf-SQUID model with finite screening,
a symmetric-series SNAIL model in its small-geometric-loop limit, physical
flux-coil JoSIM decks and quantum internal-loss baths. Flux patterns,
unequal/internal-island SNAIL modes, large SNAIL loop inductance,
multistable branch selection, degenerate 3WM and quantum fluctuations about
a strongly depleted multitone 3WM state remain outside this implementation.

Initial exclusions: full electromagnetic layout extraction, arbitrary
multi-pump Floquet networks, strong quantum non-Gaussian effects and claims of
universal junction-process accuracy. The current circuit validation uses
synthetic benchmark parameters and is not a fabrication-process calibration.
