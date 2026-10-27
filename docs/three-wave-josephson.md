# Flux-biased rf-SQUID and SNAIL traveling-wave amplifiers

Start with the [rf-SQUID](../HOWTO.md#rf-squid-3wm) and
[SNAIL](../HOWTO.md#snail-3wm) examples. Key physical references are
[Zorin](references.md#zorin2016), [Frattini et al.](references.md#frattini2017)
and [Roudsari et al.](references.md#roudsari2023).

`twpa.three_wave` models finite chains of identical flux-biased rf-SQUIDs
or SNAILs, ground capacitance, resistive loss and optional capacitively
coupled phase-matching resonators. It calculates passive dispersion,
multiharmonic 3WM gain, coherent pump depletion and pump-only Gaussian
quantum noise. Physical JoSIM decks contain individual RCSJ junctions,
loop inductors and flux-bias coils.

```python
from twpa import three_wave as tw

c = tw.example_chain('snail', cells=128, dispersion_engineered=True)
pump = .5*c.port_impedance*(.18*c.linear_chain.critical_current)**2
r = tw.three_wave_response(c, 1e9, pump, .45e9, pump*1e-6)
q = tw.quantum_noise(c, 1e9, pump, .45e9,
                    junction_temperature=.02, ground_temperature=.02)
print(r.gain, r.idler_conversion, q.added_noise, q.quantum_limit)
```

`example_chain` is a synthetic electrical benchmark, not a fabrication
preset or reproduction of a published device. For a new design, construct
`RFSQUID` or `SNAIL`, then pass it to `ThreeWaveChain(element, Cground, cells)`.
Both types retain the full current law; they do not infer Josephson behavior
from the phenomenological KIT nonlinearity coefficients.

## Potentials and equilibrium

Write phi_bar=Phi0/(2*pi), f=2*pi*Phi_external/Phi0, and let delta denote
the phase across the rf-SQUID's single junction or the SNAIL's small junction.
EJ=phi_bar*Ic refers to the rf-SQUID junction or one large SNAIL junction.

For an rf-SQUID with beta=Lloop*Ic/phi_bar,

    U/EJ = (delta-f)^2/(2*beta) - cos(delta),
    I/Ic = sin(delta) + (delta-f)/beta.

The supported range beta<1 has a unique stable equilibrium solving
delta0+beta*sin(delta0)=f. Its finite loop inductance and screening are
included in both the model and JoSIM circuit. This is a one-junction
rf-SQUID, following the circuit class discussed by
[Zorin](https://arxiv.org/abs/1602.02650). An isolated zero-inductance
dc-SQUID has a shifted sinusoidal potential and no cubic potential derivative
at its zero-current minimum; it is not substituted for this 3WM element.

For a SNAIL with n identical large junctions and small/large critical-current
ratio alpha,

    U/EJ = -alpha*cos(delta) - n*cos((delta-f)/n),
    I/Ic = alpha*sin(delta) + sin((delta-f)/n).

The model selects a unique positive-curvature zero-current root with
large-junction phase in (-pi/2,pi/2). It requires n>=2 and 0<alpha<1/n,
and rejects an ambiguous root rather than silently selecting a hysteretic
branch. Junction capacitance scales with area and resistance inversely with
area: Cs=alpha*Cj, Rs=Rj/alpha. Identical series junctions then give
Ceff=(alpha+1/n)*Cj and Geff=(alpha+1/n)/Rj. All physical junctions within
a cell share its specified bath temperature. Internal ground capacitance
on the intermediate large-junction islands is omitted.

The SNAIL potential neglects geometric loop inductance. Its physical JoSIM
adapter uses a small explicit loop inductor and requires Lgeom/Lj_large<=1e-3.
This is a model-domain guard, not a universal error estimate: verify the
small-loop limit and conditioning for a new design. Large loop inductance,
unequal series junctions and internal-island modes need a fuller circuit
model. SNAIL-based 3WM devices and dispersion engineering are discussed by
[Roudsari et al.](https://arxiv.org/abs/2209.07551); the supplied example is
not a reproduction of that device's periodic loading or gain curve.

`coefficients(element)` returns c2,c3,c4 in the expansion

    U/EJ = U0 + c2*x^2/2 + c3*x^3/6 + c4*x^4/24 + ... .

c3 generates 3WM, while c4 and cascaded higher-order effects also shift phases.
Neither the model nor a small bare c4 is assumed universally Kerr free.
The linear inductance is phi_bar/(Ic*c2); `chain.linear_chain` exposes the
equivalent small-signal circuit used by the existing ABCD and loss routines.
`bloch` and `linear_response` accept a `ThreeWaveChain` directly.
Bloch wavenumber is rad/cell and attenuation is nepers/cell; divide by
`cell_length` to obtain values per metre.

## Classical harmonic balance

Node/branch AC phases are departures from equilibrium. The pump also
rectifies the asymmetric current law. Each branch's static displacement d
is solved from zero mean current,

    mean[I(delta0 + d + x_ac(t))] = 0.

The Jacobian includes the derivative of this DC solution, a local static
Schur complement. Omitting it would give incorrect nonlinear phase shifts.
`dc_phase` returns these branch displacements in radians.

Fourier frequencies use the two-phase lattice

    f(h,m) = h*fp/2 + m*(fs-fp/2),

with h+m even, positive physical frequencies, h through 2*`harmonics`,
and m from -`sidebands` to +`sidebands`. The first modes are pump (2,0),
signal (1,1) and idler (1,-1), so fi=fp-fs. Harmonic zero with nonzero m
includes generated difference-frequency modes when extra sidebands are
retained; true DC is handled separately. Defaults are harmonics=3,
sidebands=1 and quadrature=24 per torus phase. Increase each independently
to check convergence, especially with a strong probe. The supported grid
requires sidebands*|fs-fp/2|<fp/2, distinct nonoverlapping frequencies and
nondegenerate signal/idler.

`three_wave_response` first solves the pump state and then linearizes for a
weak signal. With `depleted=True` it solves all coherent tones together.
The undepleted residual reports the pump solve; its probe response is first order.
Use additional sidebands for depletion/compression. `maximum_phase` bounds
the sum of AC element-phase magnitudes (default 0.8 rad); DC curvature is
checked separately. Newton convergence does not certify dynamical stability.

Peak phasors use exp(-i*omega*t); all electrical parameters are SI, frequency
is Hz, power W. Gain is transmitted signal power / incident signal power.
Idler conversion is a photon-flux ratio, including fs/fi. The source is a
Thevenin source with available power P and peak voltage sqrt(8*R*P).
`transmitted_power`, `reflected_power` and `incident_power` include all
retained modes. In a lossless three-tone truncation the photon changes obey
-Delta_np=Delta_ns=Delta_ni; full harmonic solutions conserve total power.

## Quantum fluctuations and internal loss

`quantum_noise` uses the same self-consistent pump and Hessian, with the
photon normalization and independent resistor baths described in the
[quantum-noise derivation](jtwpa-quantum-noise.md). Junction, ground and
resonator loss contribute both damping and vacuum/thermal fluctuations.
Temperatures may be scalar or spatial arrays using the same conventions as
the 4WM API. SNAIL's effective junction bath assumes a common temperature
for all physical junctions in a cell.

The returned Floquet modes are the signal/idler ladders fs+k*fp and fi+k*fp
at positive frequencies. These form a closed small-signal sector about a
monochromatically driven pump and its harmonics. Pump-frequency/DC sectors
are excluded. This API does not linearize about a coherently depleted
multitone signal state or model technical noise in the ideal flux source.

Here the first quantum frequency is **signal**, followed by **idler**.
Port channels are frequency-major left/right, so the primary incoming signal
is channel 0, transmitted signal is channel 1 and outgoing idler channel 3.
This differs from the 4WM result, which also retains pump-frequency modes.
The result's `node_phase` covers only the selected quantum frequencies;
obtain the coherent pump fields from `three_wave_response`.

The result exposes annihilation/creation scattering, bath transfer matrices,
normal/anomalous and quadrature covariance, bath-family contributions,
added noise, quantum limit and diagnostics. Vacuum variance is 1/2;
input-referred added noise is C_signal/G-(n_input+1/2). The finite-gain
phase-preserving lower bound is |1-1/G|/2, approaching 1/2 for high gain.
No clipping or transfer renormalization forces compliance. Diagnostics report
both commutator defects and the smallest covariance-uncertainty eigenvalue;
the validation tools assert the corresponding quantum constraints. Canonicality alone does
not establish harmonic convergence, model-domain validity or stability.

## JoSIM and validation

`chain_netlist` generates individual physical RCSJ devices and flux loops.
Series bias coils carry a common prescribed current, with coupling to each
cell's loop. Ideal bias current has zero AC fluctuation, so the coil chain
does not mediate additional AC cell coupling. A smooth quintic PWL bias ramp
reduces startup excitation; `bias_shape='linear'` retains a linear ramp.
Capacitances are included once, in the JJ models. End shunts/resonators have
half interior admittance, matching the finite model.

```python
from twpa.backends import JoSIM
from twpa.jpa_circuit import Tone, Transient

runner = JoSIM('/path/to/josim-cli', wsl_distribution='Ubuntu')
deck = tw.chain_netlist(c, [Tone(1e9, pump, ramp_time=30e-9)],
                        Transient(1e-12, 600e-9, 400e-9, 10e-12))
run = runner.run(deck, 'validation/runs/new-pump')
```

Omit WSL for native execution. Output directories must be fresh. The
[validation report](three-wave-validation.md) records actual circuit
comparisons, quantum limits, timestep/window checks and geometric-loop
limits. JoSIM provides a classical RCSJ response reference, not quantum
vacuum noise. These Gaussian models omit quantum pump depletion,
non-Gaussian dynamics, phase slips, hysteretic branch switching and
fabrication/layout extraction.
