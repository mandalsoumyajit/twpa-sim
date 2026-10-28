# Quantum noise in 4WM Josephson chains

A classical gain calculation follows the average pump and signal waves.
A quantum-noise calculation asks how small fluctuations around those waves
reach the output. Even with no thermal noise, an unused idler port and a
cold internal loss resistor carry vacuum fluctuations. Dissipation reduces
coherent amplitudes, while its environment supplies additional fluctuations.

This implementation first solves the pumped finite circuit used for classical
gain, then linearizes its response to fluctuations. It treats them as Gaussian:
means and second moments describe the state within this approximation. The
outputs are frequency-resolved noise and correlation matrices, together with
input-referred added noise. It does not evolve an arbitrary many-body
wavefunction or infer vacuum noise from random classical JoSIM voltages.

A phase-preserving amplifier must amplify both phase components while its
output still obeys quantum commutation relations. The idler and loss
environments supply additional degrees of freedom in that transformation.
The [Caves references](references.md#caves1982) establish the physical limits;
the circuit implementation below states its normalization and checks.

For a runnable starting point see the [4WM HOWTO](../HOWTO.md#4wm-jtwpa).
The annotated [Caves](references.md#caves1982),
[Caves et al.](references.md#caves2012) and [Peng et al.](references.md#peng2022)
references explain the physical context of the noise limits and Floquet channels.

`twpa.jtwpa.quantum_noise` (also `twpa.jtwpa_noise.four_wave_noise`) calculates
Gaussian quantum fluctuations about the finite circuit's coherent operating
point. It supports uniform and resonator-loaded chains, both ports, all
retained harmonics/sidebands, junction resistance, ground resistance and
resonator loss. The default operating point is pump only.

```python
from twpa import jtwpa

chain = jtwpa.Chain(critical_current=1.3e-7, ground_capacitance=1e-12,
                junction_capacitance=1e-12, cells=256,
                junction_resistance=20000., ground_resistance=100000.)
# Internal resistors retain their chosen damping; 0.05 K sets their bath
# occupations. The signal/idler are 0.9/1.1 GHz for this 1 GHz 4WM pump.
noise = jtwpa.quantum_noise(chain, pump_frequency=1e9, pump_power=7e-14,
                        signal_frequency=.9e9,
                        junction_temperature=.05, ground_temperature=.05)
print(noise.gain_db, noise.added_noise, noise.quantum_limit, noise.excess_noise)
print(noise.commutator_error, noise.minimum_uncertainty_eigenvalue)
```

Frequencies are Hz, powers W and temperatures K. Noise is a spectral density
in photon units per unit bandwidth, not an integrated noise power. Vacuum
has symmetrized noise 1/2. `output_photons` is the normally ordered spectral
occupation, with the coherent pump/probe removed. Gain is dimensionless.
The covariance describes equal, narrow, nonoverlapping frequency bins at
the retained carrier frequencies. Sweep `signal_frequency` to sample the
noise spectrum; these bins do not represent a frequency integral.

## Circuit quantization and internal baths

### Symbols and the physical meaning of a channel

A channel is one incoming or outgoing frequency at one port. A bath channel
is a fluctuation input associated with an internal resistor. The retained
frequencies include signal, idler and any additional mixing bands selected
by the numerical truncation.

| Symbol | Meaning on this page |
| --- | --- |
| `a`, `b`, `c` | External incoming, external outgoing, and internal-bath incoming field-operator vectors |
| `a_dagger` | Creation-operator vector, the Hermitian conjugate of annihilation operators; not another coherent tone |
| `U`, `V` | External annihilation/creation transfer matrices |
| `Ub`, `Vb` | Corresponding internal-bath transfer matrices |
| `R`, `f`, `h` | Real reference resistance (ohm), physical frequency (Hz), Planck constant (J s) |
| `Ic` | Physical junction critical current (A), used to normalize circuit currents |
| `n(f,T)`, `nu` | Thermal photon occupation at temperature `T` (K), and symmetrized variance `nu=n+1/2` |
| `C`, `M` | Ordinary symmetrized and anomalous output covariance matrices; `C` is not capacitance here |
| `G`, `N_added` | Signal power gain and input-referred added noise (quanta) |

Kirchhoff's current law (KCL) requires currents to balance at each node. Its
Jacobian is the matrix of derivatives with respect to solved amplitudes.
Linearizing around the pump gives the small-fluctuation response. The periodic
pump can convert an input fluctuation into other frequency channels; this is
the Floquet mixing explained in the [numerical guide](numerical-methods.md#what-the-noise-diagnostics-check).

The classical pump and its harmonics solve the same nonlinear KCL as
`four_wave_response`. Linearizing its real Jacobian gives a real-linear map
between input and output peak phasors. Rescale each voltage wave by
sqrt(2*R*h*f) to obtain photon-flux operators:

    b = U a + V a_dagger + Ub c + Vb c_dagger.

Outgoing fluctuations contain transformed external inputs plus transformed
loss-bath inputs. Creation terms account for pair generation. Keeping only
`U` would miss that process; omitting `Ub`/`Vb` while retaining resistive
damping would leave an incomplete quantum model.

Here a contains both external ports at every retained frequency; c contains
independent equilibrium bath inputs for every internal resistor at every
retained frequency. A Norton current injection associated with conductance G
has coefficient sqrt(8*h*f*G). It enters KCL divided by Ic; the classical
Jacobian already contains the resistor's damping. This normalization follows
from a resistor modeled as a matched semi-infinite transmission-line port.
It includes vacuum noise at zero temperature as well as thermal noise.

In this current-source normalization paragraph only, `G` denotes electrical
conductance, in siemens, rather than signal power gain. A Norton source is
a current source in parallel with a resistor. It represents the resistor's
environment, with a strength fixed by the conductance that causes damping.

Junction-resistor current sources enter adjacent nodes with opposite signs.
Ground-resistor conductance has half weight at either end, matching the
finite circuit. A resonator's internal resistor is eliminated together with
its node, giving an equivalent bath at the chain node with conductance

    G_eff(f) = |Yc/(Yc + Yparallel)|^2 / R_resonator.

`Yc` is the coupling capacitor's admittance, `Yparallel` is the internal
parallel resonator's admittance, and `R_resonator` is its loss resistance.
`G_eff` is the conductance seen at the chain node, in siemens. The ratio in
absolute-value bars describes how the coupling network filters that
resistor's fluctuation current before it reaches the node.

This equals the dissipative part of the effective resonator admittance.
The end resonators have half this conductance. Each eliminated bath's complex
filter phase can be absorbed into its independent equilibrium input, without
changing its covariance or commutator. Ground and resonator baths remain
separate even though both enter node KCL.

An adjoint solve for selected output rows avoids a dense circuit inverse
and one linear solve per internal resistor. The returned matrices retain all
bath channels, so memory still grows with cell count and retained modes.

`port_temperature` broadcasts to `(number_of_frequencies, 2)`: a scalar
sets all ports, a pair sets left/right, or a full array sets individual
frequency/port inputs. `junction_temperature` accepts a scalar or one value
per junction. `ground_temperature` and `resonator_temperature` accept scalars
or one value per node, including the end nodes. Temperatures set occupations
n(f,T)=1/[exp(h*f/(kB*T))-1]; they do not change the circuit parameters.
Here `kB` is Boltzmann's constant, `T` is temperature in kelvin, and the
zero-temperature limit gives `n=0`. The symmetrized variance remains `1/2`;
zero thermal occupation does not mean zero vacuum noise.

## Results and quantum constraints

In the matrix equations below, `dagger` means conjugate transpose,
`transpose` means ordinary transpose, and angle brackets mean expectation
value. Curly braces denote an anticommutator: `{A,B}=AB+BA`. `nu_port` and
`nu_bath` are diagonal matrices of each independent input's thermal-plus-
vacuum variance. Their entries are `n+1/2`, not occupation `n` alone. The
symbol `M` here is a correlation matrix, not a mutual inductance.

Port channels are frequency-major: left/right pump, left/right signal,
left/right idler, followed by extra retained modes. Thus the transmitted
primary signal is output 3 and its left input is channel 2. `frequency` and
`mode_indices` identify the retained bands. `bath_labels` identifies each
column of `bath_annihilation` and `bath_creation` as
`(family, frequency_index, junction_or_node_index)`.

The result exposes U/V (`annihilation`, `creation`), Ub/Vb,
`covariance` C=<{b_i,b_j_dagger}>/2, and anomalous covariance
M=<{b_i,b_j}>/2. For independent thermal inputs with nu=n+1/2,

    C = U nu_port U_dagger + V nu_port V_dagger
        + Ub nu_bath Ub_dagger + Vb nu_bath Vb_dagger,
    M = U nu_port V_transpose + V nu_port U_transpose
        + Ub nu_bath Vb_transpose + Vb nu_bath Ub_transpose.

`quadrature_covariance` uses `[all x, all p]`, with [x,p]=i and vacuum I/2.
For an output operator `b`, the field quadratures are
`x=(b+b_dagger)/sqrt(2)` and `p=(b-b_dagger)/(i*sqrt(2))`.
The commutator `[x,p]=xp-px=i` fixes normalization; `i=sqrt(-1)`.
`I` is an identity matrix, so vacuum covariance `I/2` means one half in each
independent quadrature. Anomalous covariance carries phase-sensitive pair
correlations between two annihilation operators.
`covariance_contributions` separates the port, junction, ground and resonator
contributions to C. These are output noise contributions, not individually
input-referred added-noise figures.

For a phase-preserving signal response, G=|U[3,2]|^2 and

    N_added = C[3,3]/G - (n_input_signal + 1/2).

`C[3,3]` is the transmitted signal's output noise variance in this ordering.
`n_input_signal` is its incoming thermal occupation. Dividing by `G` asks
what input noise would produce that output through a noiseless gain element;
subtracting the actual incoming noise leaves the added contribution.

The finite-gain vacuum lower bound is |1-1/G|/2, which becomes
(1-1/G)/2 for an amplifier and approaches 1/2 only at high gain. For a passive
attenuator with power transmission eta it is (1-eta)/(2*eta). This convention
is consistent with the [phase-preserving amplifier limit of Caves et al.](https://arxiv.org/abs/1208.5174).
`quantum_limit` reports this bound and `excess_noise` reports N_added minus
the bound. Reflection, unused conversion channels, internal loss and thermal
baths can produce excess noise even when the underlying map is canonical.

The code independently checks both canonical identities:

    U U_dagger - V V_dagger + Ub Ub_dagger - Vb Vb_dagger = I,
    U V_transpose - V U_transpose + Ub Vb_transpose - Vb Ub_transpose = 0.

It reports their largest absolute defects and the smallest eigenvalue of
the uncertainty matrix. Here `Omega` is the canonical quadrature commutator
matrix `[[0,I],[-I,0]]` for the `[all x,all p]` ordering. Physical covariance
requires the Hermitian matrix below to have nonnegative eigenvalues:
`quadrature_covariance + i*Omega/2`. A commutator defect above
`commutator_tolerance` (default 1e-6) raises an error. No matrix or covariance
is renormalized or clipped to force a quantum-limit result. Numerical
canonicality is necessary but does not establish harmonic convergence or
dynamical stability. Additional Floquet channels can change predicted noise;
see the discussion of loss ports and mode leakage in
[Peng et al.](https://doi.org/10.1103/PRXQuantum.3.020306).

## Depleted coherent states and scope

Use `depleted=True, signal_power=...` to first solve a coherently depleted
multitone state and then calculate its Gaussian differential response.
`node_phase` and `residual` describe that operating point. A strong coherent
signal generally makes its differential response phase sensitive: U and V
can both couple that same input signal to the output. In that case
`phase_preserving` is false, and `added_noise`, `quantum_limit` and
`excess_noise` are NaN. Use the returned scattering and quadrature covariance.
The gain field remains |U[3,2]|^2; it is not the strong-tone coherent power
ratio or the gain of either particular quadrature.

These calculations linearize the quantum fluctuations. They do not solve
non-Gaussian quantum dynamics, spontaneous-emission back-action on the
coherent pump, phase slips, or quantum-induced switching. The classical
solution is not dynamically certified stable. Use stable operating regimes,
check harmonic/sideband/quadrature convergence, and keep noise-induced pump
depletion negligible. The model includes the resistive losses of `Chain`;
it does not import the separate KIT/TLS noise model into a JJ circuit.

## Validation

`tests/test_jtwpa_noise.py` checks ideal finite/high-gain limits; passive
lossless scattering; cold attenuation; all three loss families; finite
temperature equilibrium and resistor emission; classical signal/idler
response; canonical commutators; covariance uncertainty; thermal spatial
profiles; coherently depleted phase-sensitive response; and spectral
resolution. Omitting internal bath channels is explicitly shown to violate
the output commutator in a lossy circuit.

Reproduce the larger circuit benchmarks with:

```powershell
python tools/validate_jtwpa_noise.py --output validation/runs/quantum.json `
  --plot validation/runs/quantum.png
```

The [quantum validation results](jtwpa-quantum-validation.md) give actual
numbers and settings. These checks use quantum identities and analytic
limits. The existing JoSIM comparisons validate the classical circuit
response underlying this Jacobian; classical JoSIM noise is not a vacuum
noise reference or an independent quantum validation.
