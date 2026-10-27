# Resonant JPA models

`twpa.jpa` implements a single-mode resonator pumped through the signal port
(Kerr JPA) and an effective flux-pumped oscillator. Both expose complex
signal/idler scattering, internal-loss noise and stability. A classical
pump-plus-probe transient solver supports finite-signal response studies.

## Units and conventions

* Resonance, pump and signal frequencies: Hz.
* `kappa_external`, `kappa_internal`, `kerr`, `pair_drive`: angular rates,
  conventionally rad/s. Convert a linewidth in Hz using `2*pi`.
* `abs(amplitude)**2`: intracavity photons.
* Incident/outgoing envelopes: square root of photons/second.
* Incident powers: watts, with photon flux `P/(h*f)`.
* Rotating frames use `exp(-i*omega*t)`, with the input-output convention
  `a_out = sqrt(kappa_external)*a - a_in`.

This frame convention is local to JPA input-output modeling. The passive TWPA
ABCD code retains its documented `exp(+i*omega*t)` convention. Their complex
phases cannot be directly interconnected without a convention/reference-plane
adapter; this release does not implement one.

## Kerr pump operating points

The model evolves a driven nonlinear resonator. Its parameters are effective
resonance frequency, total linewidth and Kerr coefficient, not junction Ic and
layout dimensions. `pump_points` finds all nonnegative roots of the normalized
pump cubic and reconstructs their complex amplitudes. Each root carries
`stable`, `unstable`, or `marginal` status from the linearized drift matrix.
Marginal means the eigenvalue real part is within `1e-9*kappa` of zero.

```python
from twpa import jpa
import numpy as np

r = jpa.Resonator(frequency=5e9, kappa_external=2*np.pi*10e6,
                  kappa_internal=2*np.pi*.25e6, kerr=-2*np.pi*10e3)
points = jpa.pump_points(r, frequency=4.994e9, power=2e-13)
for p in points:
    print(p.photons, p.stability)
p = jpa.select_pump_point(points, branch='lower')
response = jpa.kerr_response(r, p, p.frequency+np.array([-.5e6, .5e6]))
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
response = jpa.flux_response(r, pump_frequency=10e9,
                             pair_drive=.3*r.kappa, signal_frequency=[5e9])
q = jpa.quadratures(response)
print(q.power_gains, q.noise_variances)
```

Here the rotating frame is **half** the physical pump frequency. `pair_drive`
is the complex coefficient lambda in `da/dt = ... + lambda*a^dagger`.
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
trace = jpa.kerr_transient(r, p, signal_frequency=p.frequency+.5e6,
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
search, and JoSIM cross-validation remain future work.

## Example and current limits

Run `python examples/resonant_jpa.py` for an illustrative gain/noise summary.
Add `--plot gain.png` after installing `.[plot]` to save a spectrum figure.

This release adds an effective single-mode model, not a full SQUID circuit
extractor. Junction capacitance, loop screening, multimode effects and mapping
physical flux drive into the effective parameters require circuit models.
The planned JoSIM backend has not been implemented or used to validate these
JPA results. The existing traveling-wave solvers remain separate APIs.
