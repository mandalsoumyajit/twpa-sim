# Units, conventions and model boundaries

This page is a normalization reference for comparing calculations. For the
physical devices and a first worked path, begin with the
[overview](physical-overview.md) and [HOWTO](../HOWTO.md). In this page `f_p`,
`f_s` and `f_i` mean pump, signal and idler frequencies, also written `fp`,
`fs` and `fi` in code-style formulas. `omega=2*pi*f` is angular frequency;
`t` is time and `i=sqrt(-1)` sets the phasor convention.

`S11` denotes reflected wave amplitude divided by incoming amplitude at
port 1, and `S21` denotes transmitted amplitude at port 2 relative to that
input. For the stated real reference impedances, squared magnitudes give
power ratios. A Bloch wavenumber is phase advance per repeated cell; a
linewidth is the frequency width of a resonance. The JPA decay rate `kappa`
sets the linear power linewidth `kappa/(2*pi)` in Hz.

All APIs use SI electrical units. Frequency inputs are Hz; powers are watts;
temperature is kelvin. A linewidth or coupling **rate** may instead be angular:
the resonant JPA API explicitly uses rad/s for kappa, Kerr and pair drive.

## Before comparing results

| Quantity | KIT propagation | Resonant JPA | Finite Josephson chains |
| --- | --- | --- | --- |
| Complex time convention | Passive ABCD uses exp(+i*omega*t) | exp(-i*omega*t) | exp(-i*omega*t) |
| Nonlinear drive | Peak current in A | Available incident power in W / effective pair rate | Available incident power in W |
| Spatial coordinate | Cell index; extended-zone k in rad/cell | Single mode | Finite nodes; principal-zone k in rad/cell |
| Capacitance | `Cell.capacitance` per finger; two per cell | Effective mode C or circuit-adapter values | Cg per interior node, Cj per junction/effective element |
| End/reference behavior | Gain coupled-mode equations are forward only | Reflection input-output relation | Resistive two-port boundaries with half end shunts |
| Gain/noise arrays | Signal-frequency arrays | Signal-frequency arrays | One signal point per solve; scalar summary, multimode matrices |

Do not directly interconnect complex scattering matrices across time
conventions/reference planes. Conjugation and compatible port definitions
are needed; the package does not supply a universal model-interconnection adapter.

For Josephson chains, `cell_length` converts rad/cell or nepers/cell to values
per metre. It does not add an electromagnetic layout model. For KIT designs,
sections/repeats define an electrical loading pattern; a folded routing name
does not automatically insert bends, tapers or fabricated geometry.

## Frequency relations

| Process | Idler / frame |
| --- | --- |
| KIT and rf-SQUID/SNAIL 3WM | fi=fp-fs |
| Junction-chain 4WM and Kerr JPA | fi=2*fp-fs |
| Effective flux JPA | Frame=physical pump/2; fi=2*frame-fs |

Require positive supported sideband frequencies. Distributed Josephson solvers
currently require distinct signal/idler bands. For degenerate resonant JPA
operation, use quadratures; treating one mode as two independent vacuum inputs
would give the wrong noise prediction.

## Power and photon flux

For a real port impedance R and peak voltage/current waves:

`R` is measured in ohms; `V_wave` and `I_wave` are complex **peak**, rather
than root-mean-square, voltage and current amplitudes. `P` is average wave
power in watts, `h` is Planck's constant and `f` is physical frequency in Hz.
The unit dBm measures power relative to 1 milliwatt.

    P = |V_wave|^2/(2*R) = R*|I_wave|^2/2
    photon_flux = P/(h*f)
    P[W] = 1e-3 * 10**(P[dBm]/10)

The Thevenin source peak voltage for available incident power P is
sqrt(8*R*P); the incoming wave is half that voltage. Actual branch currents
depend on the solved circuit and reflections. Example powers constructed
from a reference current do not impose that current in every junction.

Signal power gain equals signal photon-flux gain because input/output signal
frequencies coincide. Josephson `idler_conversion` is a **photon-flux** ratio:
multiply by fi/fs for an idler power ratio. KIT nonlinear `currents` are
peak amperes in the order pump, signal, idler, rather than photon amplitudes.

In 3WM chains, `linear_chain.critical_current` is Ic*c2, a curvature-derived
linear scale. The physical element's `critical_current` is the single
rf-SQUID or large SNAIL junction Ic. They are not interchangeable.

## Quantum noise

Here `G` means signal power gain, `n` means thermal photon occupation, and
`I` in matrix expressions is an identity matrix rather than electrical
current. A dagger denotes a creation operator or conjugate transpose as
appropriate to the object. A covariance contains second moments of field
fluctuations; it excludes the coherent pump or signal amplitude.

Symmetrized spectral vacuum noise is 1/2 quantum. Gaussian quadratures use
[x,p]=i and vacuum covariance I/2. Normally ordered occupation subtracts that
half quantum; coherent tones are excluded from fluctuation covariance.

Input-referred added noise subtracts the incoming signal noise after dividing
output signal noise by gain. The phase-preserving vacuum lower bound is
|1-1/G|/2: (1-1/G)/2 for G>1, and (1/G-1)/2 for attenuation. Half a quantum
is the high-gain limit, not the exact finite-gain bound.

Internal resistance causes damping **and** bath fluctuations. Omitting its
vacuum input can violate output commutators. Temperature changes bath
occupation without automatically changing Ic, R or C. The KIT/TLS and
Josephson noise pipelines are separate models, not one universal loss adapter.

| Quantum result | First modes / scalar noise scope |
| --- | --- |
| `jpa` | Nambu vector [signal, idler-dagger]; use quadratures at degeneracy |
| `jtwpa.quantum_noise` | Pump, signal, idler, then retained modes; frequency-major left/right ports |
| `three_wave.quantum_noise` | Signal, idler, then Floquet ladders; frequency-major left/right ports; pump-only operating point |

A Nambu vector groups annihilation and creation components, here signal and
conjugate idler, because the pump couples them. A Floquet ladder groups
frequencies linked by the periodic drive. These are choices of representation,
not additional physical junctions or ports.

For 4WM the primary outgoing signal channel is 3 and incoming left signal is
2. For 3WM they are 1 and 0. Prefer named result summaries unless you need
the full matrices. See [quantum noise](jtwpa-quantum-noise.md) and
[3WM noise](three-wave-josephson.md#quantum-fluctuations-and-internal-loss).

## Model boundaries

- KIT `coupled`: forward small-signal gain/noise with pump attenuation;
  conversion-induced depletion uses the separate classical `nonlinear` API.
- JPA: effective single mode; pump roots have explicit linear stability.
  A physical flux-to-pair-rate conversion needs a circuit adapter.
- 4WM: finite uniform/RPM RCSJ chain; Gaussian noise about a classical
  operating point. A coherent depleted state can be phase sensitive.
- 3WM: unique-branch rf-SQUID or identical-series SNAIL in the small-loop
  limit; coherent depletion is available, quantum noise is pump-only.

Distributed harmonic-balance convergence does not certify dynamical stability.
These propagation/noise models do not perform layout extraction, quantum-induced switching or
general non-Gaussian dynamics. Consult the specific guide before extending a regime.

The separate [layout-to-circuit workflow](layout-backannotation.md) assembles
finite Josephson circuits from supported field extraction, and KIT periodic
circuits from explicitly supplied calibration. These are distinct APIs with
recorded parameter origins; extraction does not automatically modify every
propagation or quantum API listed above.

## Troubleshooting

| Symptom | Check |
| --- | --- |
| No selected stopband / no pump candidates | Use an explicit pump for an unloaded line; check grid range and all sidebands |
| Pump or sideband outside passband | Inspect passive dispersion before changing nonlinear controls |
| Effective-stub pole on a grid point | Resolve and avoid the pole itself; do not smooth a singularity into a passband |
| Multiple JPA roots | Choose a stable branch or use documented continuation |
| NaN resonant added noise at degeneracy | Use `jpa.quadratures` |
| JJ Newton failure / phase guard | Reduce drive; refine continuation and spectral resolution before interpreting gain |
| Commutator defect | Check pump convergence/quadrature; retain physical loss baths; do not renormalize away the defect |
| JoSIM failure | Inspect logs, executable/capabilities and the fresh output path |
| JoSIM gain drifts with window | Extend turn-on/settling/acquisition; resolve the highest relevant oscillations |

The [HOWTO](../HOWTO.md) shows the complete workflow. Physical background
and software references are collected in [References](references.md).
