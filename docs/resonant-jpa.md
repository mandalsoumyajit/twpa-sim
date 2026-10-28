# Resonant JPA models

A resonant Josephson parametric amplifier stores microwave energy in an
inductor–capacitor mode whose inductance is nonlinear. A strong pump changes
that mode's response so a weak incident signal can emerge with larger power.
The examples here are reflection amplifiers: the outgoing signal is measured
at the same external port used to drive the resonator. The pump supplies
energy to the signal and its partner frequency, the idler.

There are several useful starting levels. `twpa.jpa` takes an effective
resonance, decay rates and nonlinear coefficient for fast operating-point
and noise studies. The physical `ParallelJPA`/`SQUIDJPA` adapters start from
junctions, capacitors, inductors and a source termination. The
[layout extraction](magnetic-extraction.md) derives supported capacitance and
interconnect parameters from drawn geometry. These are different levels of
description, not interchangeable parameter lists.

Start with the [Kerr](../HOWTO.md#kerr-jpa) and
[flux-pump](../HOWTO.md#flux-pumped-jpa) HOWTO examples. Related literature:
[Eichler–Wallraff](references.md#eichler2014) and [Caves](references.md#caves1982).

`twpa.jpa` implements a single-mode resonator pumped through the signal port
(Kerr JPA) and an effective flux-pumped oscillator. Both expose complex
signal/idler scattering, internal-loss noise and stability. A classical
pump-plus-probe transient solver supports finite-signal response studies.

## Units and conventions

| Symbol | Meaning | API field or units |
| --- | --- | --- |
| `f_0`, `f_p`, `f_s`, `f_i` | Unpumped resonance, pump, signal and idler frequencies | Hz; `frequency` sets `f_0` |
| `omega=2*pi*f` | Angular frequency corresponding to ordinary frequency `f` | rad/s |
| `a` | Complex intracavity amplitude; `abs(a)**2` is mean cavity photon number | Square root of photons |
| `a_in`, `a_out` | Incoming/outgoing amplitudes at the external port | Square root of photons/s |
| `kappa_e`, `kappa_i`, `kappa` | External, internal and total energy-decay rates; `kappa=kappa_e+kappa_i` | `kappa_external`, `kappa_internal`; rad/s |
| `K` | Angular resonance shift per intracavity photon | `kerr`; rad/s per photon |
| `Delta=2*pi*(f_p-f_0)` | Angular pump detuning; positive means pump above resonance | rad/s |
| `P_p`, `h` | Incident pump power and Planck's constant | W; J s |

The exponential time convention uses `i=sqrt(-1)` and time `t` in seconds.
It fixes complex phases, not measured power gain. External decay means
energy escaping through the measurement port; internal decay means energy
escaping into an unobserved environment. A linear mode's power linewidth in
Hz is `kappa/(2*pi)`.

Rotating frames use `exp(-i*omega*t)`, with the input-output convention
`a_out = sqrt(kappa_external)*a - a_in`. This equation combines radiation
from the stored cavity field with reflection of the incoming wave. Power
and photon flux differ by photon energy `h*f`; the API accepts incident
powers in watts and performs that conversion internally.

This frame convention is local to JPA input-output modeling. The passive TWPA
ABCD code retains its documented `exp(+i*omega*t)` convention. Their complex
phases cannot be directly interconnected without a convention/reference-plane
adapter; this release does not implement one.

## Kerr pump operating points

In the pump's rotating frame, the classical effective equation is

    da/dt = (i*Delta - kappa/2)*a - i*K*abs(a)^2*a + sqrt(kappa_e)*a_in.

The rotating frame removes the rapid oscillation at the pump frequency so
`a` varies slowly. The drive has `abs(a_in)**2=P_p/(h*f_p)` photons per second.
At a steady pump state the photon number `n=abs(a)**2` therefore obeys

    n * [(kappa/2)^2 + (Delta-K*n)^2] = kappa_e*P_p/(h*f_p).

This cubic expresses feedback between stored energy and resonance shift.
It can have several nonnegative solutions. Linear stability asks whether a
small perturbation decays; finding a root alone does not answer that question.

`pump_points` finds all nonnegative roots of the normalized
pump cubic and reconstructs their complex amplitudes. Each root carries
`stable`, `unstable`, or `marginal` status from the linearized drift matrix.
Marginal means the eigenvalue real part is within `1e-9*kappa` of zero.

```python
from twpa import jpa
import numpy as np

resonator = jpa.Resonator(frequency=5e9, kappa_external=2*np.pi*10e6,
                  kappa_internal=2*np.pi*.25e6, kerr=-2*np.pi*10e3)
points = jpa.pump_points(resonator, frequency=4.994e9, power=2e-13)
for candidate in points:
    print(candidate.photons, candidate.stability)
# Choose a stable occupation before perturbing the pump with a weak signal.
operating_point = jpa.select_pump_point(points, branch='lower')
signal_frequencies = operating_point.frequency + np.array([-.5e6, .5e6])
response = jpa.kerr_response(resonator, operating_point, signal_frequencies)
print(response.gain_db, response.added_noise)
```

In bistability, selection without a branch raises an error. `lower` and `upper`
refer to occupation, not frequency. Passing `previous=last_point` instead
chooses the nearest stable complex amplitude for discrete sweep continuation.
This does not simulate noise-induced switching or guarantee branch tracking
across coarse sweep steps. Small-signal spectra reject unstable/marginal points
and pump points created for another resonator.

The implementation follows the driven Kerr-resonator input-output treatment of
[Eichler and Wallraff, EPJ Quantum Technology 1, 2 (2014)](https://link.springer.com/article/10.1140/epjqt2),
including the steady-state cubic, fluctuation matrix and internal-loss ports.
Tests compare the resolvent with the closed scattering expression, verify
stable/unstable branches and check pump-equation residuals. They do not claim
experimental calibration.

## Gain and noise away from degeneracy

The pump state is solved first. The signal is then small enough that only
its first-order perturbation is retained. This is why the result is a
small-signal gain curve rather than a finite-input compression curve.
In the operator notation below, a dagger denotes Hermitian conjugation,
so `a_i^dagger` is an idler creation operator. The vector groups the signal
and conjugate idler because parametric amplification couples those two
fluctuations. The scattering matrix `S` maps incoming to outgoing fields;
its row/column indices start at zero in Python.

For the Kerr model, `frame_frequency` is the pump frequency. For a signal at
`f_s`, the conjugate idler is at `2*frame_frequency - f_s`. `scattering` acts
on `[a_s, a_i^dagger]`; `bath_scattering` acts on the corresponding internal
bath inputs. `gain` is `abs(S[0,0])**2`, and `conversion_gain` is
`abs(S[0,1])**2`. These are photon-flux coefficients. Signal power gain equals
its photon-flux gain; converting idler flux to power needs its frequency ratio.

`added_noise` is input-referred symmetrized quanta. Signal-input noise is
excluded, the external idler is vacuum, and internal bath occupations follow
the supplied temperature at each physical sideband frequency. With no internal
loss, gain minus conversion is one and added noise is `(1-1/G)/2`. At zero
signal gain the input-referred quantity is undefined/infinite.

## Flux pump and degenerate quadratures

```python
response = jpa.flux_response(resonator, pump_frequency=10e9,
                             pair_drive=.3*resonator.kappa, signal_frequency=[5e9])
quadratures = jpa.quadratures(response)
print(quadratures.power_gains, quadratures.noise_variances)
```

Here the rotating frame is **half** the physical pump frequency. `pair_drive`
is the complex coefficient lambda in `da/dt = ... + lambda*a^dagger`.
Here `lambda` is a pair-coupling rate, not the London penetration length
used in magnetic extraction. Its magnitude sets the effective pump strength
and its phase selects the amplified quadrature axis.
It is an input, not a conversion from pump power or applied flux. An optional
`frequency_shift` in Hz supplies an effective pump-induced resonance shift.
The unpumped zero-mean state is accepted only below its linear instability
threshold. Kerr saturation and operation above threshold are not modeled by
this effective flux-pump solver, even if the Resonator contains a nonzero Kerr.

At degeneracy, signal and idler are one mode. `added_noise` is deliberately
NaN; a phase-preserving half-photon formula would be inappropriate.
`quadratures` requires one degenerate frequency point and returns real
quadrature transfer matrices, symmetrized output covariance for vacuum external
input, singular-value power gains, and covariance eigenvalues. Vacuum variance
is 1/2. These are **output noise spectra at zero analysis frequency**, not
intracavity steady-state variances. Nonzero signal offsets require the
two-sideband scattering representation instead.

## Finite-signal transients

`kerr_transient` starts at a chosen stable pump state and integrates the
classical cavity equation with an added coherent probe. For example:

```python
trace = jpa.kerr_transient(resonator, operating_point,
                          signal_frequency=operating_point.frequency+.5e6,
                          signal_power=1e-18, duration=20e-6)
print(trace.signal_gain, trace.idler_conversion)
```

The routine discards an initial fraction of the record and least-squares fits
DC and the first three positive/negative probe harmonics. It requires at least
two retained beat periods and checks output sampling for harmonic aliasing.
The integration step also resolves the cavity lifetime and probe beat.
Vary duration, discard fraction, tolerances and sampling for a compression
study. Strong-signal harmonic truncation and slow switching near a bifurcation
need extra convergence checks. Transients are classical; they do not sample
quantum noise. Dynamic pump sweeps, a dedicated resonant 1-dB compression
search remain future work. [JoSIM circuit comparisons](josim-validation.md)
now check classical Kerr gain, idler conversion and finite-signal compression.

## Example and current limits

Run `python examples/resonant_jpa.py` for an illustrative gain/noise summary.
Add `--plot gain.png` after installing `.[plot]` to save a spectrum figure.

The effective single-mode model assumes its parameters are already known.
Deriving them from junction capacitance, loop screening or a physical flux
drive requires a circuit model, and extra modes can limit a single-mode
description. The [JoSIM backend](josim-backend.md) supplies physical RCSJ and
mutually coupled SQUID circuits. The [layout workflow](magnetic-extraction.md)
adds the dedicated single-junction JPA extraction. The
[general terminal-network adapter](layout-backannotation.md) retains both
SQUID junctions and finite loop/return inductance, with explicitly prescribed
flux and its own classical response. It does not silently fit that circuit
to the effective single-mode parameters.

The validation reports compare classical Kerr and flux-pump response,
including quadrature gains, in specified weak-phase and weak-modulation
regimes. Quantum noise is checked against analytic identities rather than
classical transients. Choose a traveling-wave guide for spatial propagation
and finite-chain effects, which require separate APIs.
