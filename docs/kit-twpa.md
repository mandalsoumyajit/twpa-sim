# KIT ladders, nonlinear propagation and TLS loss

A kinetic-inductance traveling-wave amplifier uses the current dependence of
a superconducting film's inductance to exchange energy between propagating
microwaves. A strong pump and a weak signal enter a long line; repeated cells
allow the signal to build up as it travels. With a direct-current bias, the
main model supports three-wave mixing: pump energy supplies a signal–idler
pair. This is the physical setting discussed by
[Malnou et al.](references.md#malnou2021), with earlier kinetic-inductance
amplifier background in [Eom et al.](references.md#eom2012).

The simulation has three stages. First, the electrical ladder determines
dispersion and passive reflection/transmission. Next, a propagation model
calculates pump, signal and idler amplitudes using the chosen current
nonlinearity. Finally, an optional dielectric-loss model adds pump attenuation
and noise from the lossy material. The cells are electrical equivalents;
their capacitance and inductance need calibration or extraction for a real
device. The separate layout generator does not automatically supply that
calibration.

## Quantities used in this guide

| Symbol | Physical meaning | Units |
| --- | --- | --- |
| `f_p`, `f_s`, `f_i` | Pump, signal and idler frequencies; for 3WM, `f_i=f_p-f_s` | Hz |
| `I_dc`, `I_scale` | Applied direct current and characteristic current controlling inductance nonlinearity | A |
| `epsilon`, `xi` | Effective three-wave and Kerr coefficients returned by the bias law | 1/A; 1/A² |
| `kerr_factor` | Convention/override multiplier in the effective Kerr coefficient | Dimensionless |
| `k` | Phase advance per cell from Bloch analysis of the periodic ladder | rad/cell |
| `S11`, `S21` | Complex reflected and transmitted wave amplitudes relative to the incident wave | Dimensionless |
| `G` | Forward signal power gain | Dimensionless; `10*log10(G)` is gain in dB |

An electrical cell has series inductance and two capacitive fingers. A Bloch
calculation asks how a wave changes from one repeated cell to the next. A
finite scattering calculation includes the actual line length and terminations.
See [the numerical introduction](numerical-methods.md#passive-circuits-before-pumping)
for the distinction.

The top-level `twpa` APIs model periodic electrical ladders and a
kinetic-inductance current law, with an optional phenomenological 3WM/Kerr
coefficient pair. Start with [HOWTO examples 01 and 02](../HOWTO.md#kit-and-tls-loss).
Background is given by [Eom et al.](references.md#eom2012) and
[Malnou et al.](references.md#malnou2021); the included circuits are not
claimed to reproduce their experimental devices.

## Build a design

`Cell` specifies inductance, capacitance per finger, optional stub inductance,
dielectric participation and topology. There are two capacitive fingers per
cell. `Section(cell, count)` groups identical cells; `Design` repeats a tuple
of sections and optionally appends a tail. Values are electrical equivalents,
not an electromagnetic extraction of a layout.

`preset(name)` returns a design and bias. `v1-target`, `v1-cpw`, `v1-idc` and
`v2` are historical electrical alternatives. `v2-tls`/`v2-standard`,
`v2-folded` and `v2-four-pass` share a calibrated cell model with different
lengths. They do not automatically simulate bends or four-pass routing.
See the [source audit](source-audit.md) for parameter and Kerr-factor history.

`Bias` takes peak pump/DC currents, an inductance scaling current and port
impedance. Without an override the current law supplies

    epsilon = 2*I_dc/(I_scale^2 + I_dc^2)
    xi = kerr_factor/(I_scale^2 + I_dc^2).

An explicit `Nonlinearity(three_wave=..., kerr=...)` supplies coefficients in
1/A and 1/A^2. This changes the envelope interaction; it does not add JJ
plasma response or a microscopic Josephson potential.

## Passive and nonlinear calculations

`linear(design, frequency, impedance)` returns finite S11/S21, extended-zone
Bloch wavenumber, attenuation, passbands and stopband intervals. Wavenumber
is rad/cell and attenuation nepers/cell. Resolve narrow loading features and
avoid sampling an effective-stub pole exactly. `network_scattering` cascades
arbitrary electrical blocks without assigning them one homogeneous gain law.

`small_signal` uses a forward undepleted 3WM solution with Kerr phase shifts.
Provide pump frequency explicitly for a uniform line; automatic selection
searches near a selected loading stopband. For 3WM, `f_i=f_p-f_s`. Dispersion
and nonlinear phase shifts must also let contributions from successive cells
add coherently; energy conservation alone does not establish phase matching.

`nonlinear` integrates pump, signal and idler peak currents together. It
supports fractional inductance gradients and normalized routing profiles;
`route_profile` provides simple examples, not a full layout solver.
`compression` sweeps input power and finds the first bracketed 1 dB crossing,
returning NaN if none was bracketed. These finite-signal routines are classical.

## TLS loss and noise

`TLS` supplies an intrinsic low-power loss tangent, temperature, critical
field, field/current conversion and linewidth. `tls_tangent` combines thermal
tanh suppression, saturation and a phenomenological Lorentzian pump hole.
The hole prevents a distant pump from unrealistically saturating every TLS
frequency. Its linewidth and material parameters must be calibrated.
General dielectric/TLS motivation is discussed by
[Martinis et al.](references.md#martinis2005) and
[Gao et al.](references.md#gao2008); these references do not prescribe this
repository's device-specific participation or saturation-hole parameters.

`idc_fields` solves a finite-volume electrostatic half-gap edge domain.
`oxide_design` compares coated and bare field energies, rescales calibrated
cell capacitances by their ratio, and assigns oxide participation. Energy
ratios are useful here because this small domain is not an extraction of
an entire IDC or package. Check mesh, air/substrate domain and oxide-thickness
resolution before drawing quantitative conclusions.

`coupled` feeds local lossy Bloch parameters into evolution of pump current,
signal/idler transfer and thermal bath covariance. Pump attenuation is
included by default; conversion-induced pump depletion is not. Gain is
forward only, with no reflected-wave feedback into the amplifier's
coupled-mode propagation equations.
The returned noise is input referred, symmetrized and photon normalized,
with a vacuum input idler. `constant_langevin` provides the constant-parameter
matrix-exponential solution and bath integral for comparison.
The loss is represented by effective thermal baths, not a microscopic
nonequilibrium driven-TLS noise calculation.

Participation weights how strongly the electric field samples the lossy
oxide. A loss tangent describes dielectric dissipation; the critical field
sets the scale at which this phenomenological TLS loss saturates.
`field_per_amp` converts line current to an estimated local electric field.
These quantities need material/geometry calibration independently of the
superconducting film's nonlinear current scale.

## Reading a gain and loss calculation

In [example 01](../examples/howto/01_kit.py), compare undepleted weak-signal
gain with finite-signal gain. The latter permits the signal to draw energy
from the pump. In [example 02](../examples/howto/02_tls_loss.py), inspect the
electric-energy participation and capacitance ratio before looking at gain.
The coupled solver then reports what remains of the pump and the bath's
contribution to input-referred added noise. A small electrostatic domain can
converge numerically without representing the full capacitor or package;
keep that model boundary separate from mesh accuracy.

The quantum bound and [noise conventions](conventions.md#quantum-noise)
follow [Caves](references.md#caves1982). A cold lossy line must inject its
bath's vacuum fluctuations; a lossy classical transfer alone is insufficient.
The coupled pipeline remains available alongside the newer Josephson APIs.

## Verification and limitations

Tests compare source implementations, analytic ladders, passive attenuation,
lossless gain/noise limits, mesh-related pipeline behavior and weak/finite
probe consistency. The [source audit](source-audit.md) explains which original
plots used supplied tables rather than a fully computed pipeline. The current
JoSIM reports cover Josephson circuits, not this KIT/TLS implementation.

This model homogenizes periodic loading for gain and uses electrical
parameters supplied by the user. Film dispersion, material calibration,
geometry extraction, nonlinear coefficient validity and numerical convergence
remain design-specific. See [conventions](conventions.md) and
[References](references.md#kit-and-dielectric-loss).
