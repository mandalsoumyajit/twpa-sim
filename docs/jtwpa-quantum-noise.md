# Quantum noise in 4WM Josephson chains

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

c = jtwpa.Chain(critical_current=1.3e-7, ground_capacitance=1e-12,
                junction_capacitance=1e-12, cells=256,
                junction_resistance=20000., ground_resistance=100000.)
r = jtwpa.quantum_noise(c, pump_frequency=1e9, pump_power=7e-14,
                        signal_frequency=.9e9,
                        junction_temperature=.05, ground_temperature=.05)
print(r.gain_db, r.added_noise, r.quantum_limit, r.excess_noise)
print(r.commutator_error, r.minimum_uncertainty_eigenvalue)
```

Frequencies are Hz, powers W and temperatures K. Noise is a spectral density
in photon units per unit bandwidth, not an integrated noise power. Vacuum
has symmetrized noise 1/2. `output_photons` is the normally ordered spectral
occupation, with the coherent pump/probe removed. Gain is dimensionless.
The covariance describes equal, narrow, nonoverlapping frequency bins at
the retained carrier frequencies. Sweep `signal_frequency` to sample the
noise spectrum; these bins do not represent a frequency integral.

## Circuit quantization and internal baths

The classical pump and its harmonics solve the same nonlinear KCL as
`four_wave_response`. Linearizing its real Jacobian gives a real-linear map
between input and output peak phasors. Rescale each voltage wave by
sqrt(2*R*h*f) to obtain photon-flux operators:

    b = U a + V a_dagger + Ub c + Vb c_dagger.

Here a contains both external ports at every retained frequency; c contains
independent equilibrium bath inputs for every internal resistor at every
retained frequency. A Norton current injection associated with conductance G
has coefficient sqrt(8*h*f*G). It enters KCL divided by Ic; the classical
Jacobian already contains the resistor's damping. This normalization follows
from a resistor modeled as a matched semi-infinite transmission-line port.
It includes vacuum noise at zero temperature as well as thermal noise.

Junction-resistor current sources enter adjacent nodes with opposite signs.
Ground-resistor conductance has half weight at either end, matching the
finite circuit. A resonator's internal resistor is eliminated together with
its node, giving an equivalent bath at the chain node with conductance

    G_eff(f) = |Yc/(Yc + Yparallel)|^2 / R_resonator.

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

## Results and quantum constraints

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
`covariance_contributions` separates the port, junction, ground and resonator
contributions to C. These are output noise contributions, not individually
input-referred added-noise figures.

For a phase-preserving signal response, G=|U[3,2]|^2 and

    N_added = C[3,3]/G - (n_input_signal + 1/2).

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
