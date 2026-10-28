# Four-wave Josephson traveling-wave amplifiers

This model describes a microwave transmission line assembled from many
Josephson junctions in series, with capacitance from each intervening node
to ground. The junction inductance depends on phase, so a strong propagating
pump can transfer energy into a weaker signal and an idler. An unbiased
sinusoidal junction law supplies the four-wave-mixing interaction: two pump
photons provide the energy of a signal–idler pair.

Unlike an ideal forward-wave propagation model, this calculation keeps a
finite number of circuit nodes, their source/load resistors and the reflected
waves. It first finds the driven pump state, then either perturbs that state
with a weak signal or solves a finite signal and pump together. Optional
resonant loading changes the dispersion to improve phase matching, following
the physical approach of [O'Brien et al.](references.md#obrien2014).

For example, the main tutorial asks what leaves a 64-cell line when pumped
at 1 GHz and weakly probed at 0.9 GHz. Its principal idler is 1.1 GHz. The
reported gain is transmitted signal power at the far end divided by incident
signal power. This short synthetic example is useful for learning the circuit
and noise calculations, rather than demonstrating an optimized high-gain design.

Start with the [commented HOWTO](../HOWTO.md#4wm-jtwpa). Background references
include [O'Brien et al.](references.md#obrien2014) for resonant phase matching
and [Peng et al.](references.md#peng2022) for Floquet/noise context.

`twpa.jtwpa` models an unbiased series Josephson-junction chain with ground
capacitance, junction capacitance and resistive loss. `ResonantLoad` adds an
identical capacitively coupled parallel LC resonator at each node for resonant
phase matching (RPM). The finite circuit has resistive source/load ports and
half of a cell's shunt admittance at each end. These same physical parameters
generate full-sine RCSJ circuits for JoSIM.

```python
import numpy as np
from dataclasses import replace
from twpa import jtwpa
from twpa.jpa_circuit import PHIBAR

# PHIBAR is the flux quantum divided by 2*pi. I_c=PHIBAR/L_J chooses
# a zero-phase junction inductance L_J of 2.5 nH.
chain = jtwpa.Chain(critical_current=PHIBAR/2.5e-9,
                ground_capacitance=1e-12, junction_capacitance=1e-12,
                cells=384)
# Match the source/load to the periodic line's small-signal pump impedance.
chain = replace(chain, port_impedance=float(jtwpa.bloch(chain, [1e9]).impedance[0].real))
# P=R*I_peak^2/2 specifies an incoming wave, not a forced current at every JJ.
pump_power = .5*chain.port_impedance*(.4*chain.critical_current)**2
result = jtwpa.four_wave_response(chain, pump_frequency=1e9,
    pump_power=pump_power, signal_frequency=.9e9, signal_power=pump_power*1e-6)
print(10*np.log10(result.gain), result.idler_conversion, result.pump_transmission)
```

Run `python examples/four_wave_jtwpa.py` to compare uniform and RPM chains.
`linear_response` returns exact finite S11/S21 using scaled ABCD cascades.
`bloch` returns principal-zone wavenumber in rad/cell, attenuation in
nepers/cell, Bloch impedance and the corresponding lossless-cell passband
classification. Divide wavenumber by `cell_length` to obtain rad/m. Near poles
and band edges use finite scattering and adequately resolved frequency grids.

## Equations and normalization

| Symbol | Meaning | Units |
| --- | --- | --- |
| `Phi_left`, `Phi_right` | Node fluxes on either side of a junction; a node-flux derivative is its voltage | Wb |
| `delta` | Junction phase difference `(Phi_left-Phi_right)/phi_bar` | rad |
| `phi_bar`, `hbar`, `e` | Reduced flux quantum `hbar/(2*e)`, reduced Planck constant, positive elementary charge | Wb; J s; C |
| `I`, `I_c` / `Ic` | Junction branch current and physical critical current | A |
| `C_j` / `Cj`, `C_g` / `Cg` | Capacitance across a junction and capacitance from an interior node to ground | F |
| `R_j` / `Rj` | Junction shunt loss resistance | ohm |
| `f_p`, `f_s`, `f_i` | Pump, signal and idler frequencies, also written `fp`, `fs`, `fi` below | Hz |
| `R`, `P`, `V` | Real port resistance, available incident power, peak voltage-wave phasor | ohm; W; V |

An overdot means a derivative with respect to time: `delta_dot=d(delta)/dt`
and `delta_ddot=d²(delta)/dt²`. The imaginary unit is `i`, and
`omega=2*pi*f` is angular frequency. In the code, `PHIBAR` is `phi_bar`.

For node flux Phi and branch phase delta=(Phi_left-Phi_right)/phi_bar,
phi_bar=hbar/(2e), the branch current is

    I = Ic*sin(delta) + Cj*phi_bar*delta_ddot + phi_bar*delta_dot/Rj.

Ground capacitors, optional ground resistors, and optional resonator
admittances enter Kirchhoff's current law (KCL): currents into and out of each
node must balance. Admittance is current divided by voltage. The model solves this finite circuit directly in
frequency space; reflections and end loads enter the same equations as the
nonlinear branches. It does not assume a forward envelope or a continuum
dispersion relation. The series linear impedance includes the JJ plasma
response; junction capacitance is included once in both model and netlist.

Peak phasors use exp(-i*omega*t). Frequencies are Hz, powers W, all electrical
values SI. `node_phase` is Phi/phi_bar. Source power means the available
incident power into a real reference impedance R. Source peak voltage is
sqrt(8*R*P), incident voltage sqrt(2*R*P), and a port voltage phasor V carries
power |V|^2/(2*R). The opposite Fourier convention in the older KIT circuit
module must be conjugated before comparing complex phasors.

Four-wave mixing conserves fi=2*fp-fs. `gain` is transmitted signal **power**
divided by incident signal power. `idler_conversion` is transmitted idler
**photon flux** divided by incident signal photon flux; multiply it by fi/fs
for the power conversion ratio. `pump_transmission` is output pump power
divided by incident pump power. `reflection` is the complex signal reflection
coefficient. The three power arrays include every retained frequency at both
ports and permit power-balance checks.

## Harmonics, sidebands and depletion

The nonlinear sine current produces more than the three principal tones.
Harmonic balance approximates the motion by a chosen set of sinusoidal
components and solves their amplitudes self-consistently. A harmonic is an
integer multiple of the pump frequency; an offset sideband lies near such a
harmonic. Keeping too few components can change the predicted gain even
when the chosen equations have a small residual.

By default `nonlinearity='sine'` projects the exact sine current onto a
two-phase Fourier basis. Retained frequencies are

    f(h,m) = h*fp + m*(fs-fp),

where `h` and `m` are integer frequency indices: `h` labels a pump harmonic,
and `m` counts steps of the signal–pump separation `fs-fp`. The `h` in this
frequency-index formula is not Planck's constant. A two-phase Fourier basis
uses independent periodic pump and offset phases to represent these tones.

The basis keeps positive odd `h` through `harmonics` (default 3), and integer `m` from
`-sidebands` to `+sidebands` (default 1). Results start with pump, signal,
idler, followed by the extra modes; `mode_indices` gives the explicit map.
The torus quadrature size defaults to 16 per phase. Increase quadrature,
harmonics and sidebands independently until the observables converge.
The default basis includes third-harmonic mixing even for a weak signal.
At stronger probes additional offset sidebands matter.

`nonlinearity='cubic'` selects sin(delta) approximately delta-delta^3/6,
retaining only pump, signal and idler. It provides a quartic-potential 4WM
approximation with self/cross-phase shifts and pair conversion. It can be
less accurate than the full-sine multiharmonic solution well below switching.

The default undepleted calculation first solves the nonlinear pump and its
harmonics, then linearizes the real sparse Jacobian for an infinitesimal
signal. Its residual reports the pump solve; signal power scales the linear
perturbation and should remain small. `depleted=True` solves all retained
amplitudes together, including pump depletion and generated sidebands.
Sweep `signal_power` explicitly for compression. Continuation in pump drive
and Newton line search improve convergence, but do not establish dynamical
stability or select every possible branch. A conservative sum of junction
phase-tone magnitudes is returned and checked against `maximum_phase`
(default 1.2 rad). This is a truncation guard, not a switching criterion.

## Scope and validation

To interpret a run, inspect passive transmission before choosing the pump.
Then check signal gain, idler conversion and pump transmission together.
Changing cell count or pump power changes accumulated nonlinear phase as
well as interaction strength, so increased drive need not improve gain.
Resonant phase matching adjusts the passive dispersion to compensate this
phase mismatch; it does not directly impose a target gain. For stronger
signals, compare depleted and weak-signal responses and examine all retained
output powers. The [numerical-method guide](numerical-methods.md) explains why
spectral convergence, a phase guard and dynamical stability are different checks.

Uniform and identical-resonator-loaded chains, passive dispersion, gain,
idler generation, harmonics and pump depletion have been checked against
[actual JoSIM transient simulations](jtwpa-josim-validation.md). Analytic tests
also cover plasma response, finite scattering, Fourier coefficients,
Jacobian derivatives, power and three-tone Manley–Rowe balances.

`obrien_2014()` provides the circuit values from
[O'Brien et al., resonant phase matching](https://arxiv.org/html/1406.2346v1):
Ic=3.29 uA, Cg=39 fF, Cj=329 fF, coupling C=10 fF, resonator C=7.036 pF,
L=100 pH and 2000 cells of length 10 um. Junction inductance is derived from
Ic, rather than separately imposing its rounded paper value. This preset
does not claim to reproduce the paper's full continuum gain curves.

Current limits: unbiased identical junctions; distinct nearby signal/idler
bands; no phase slips, DC voltage state, arbitrary fabrication profiles,
SQUID/SNAIL flux bias, or dynamic stability calculation.
[Quantum Gaussian scattering and internal-loss noise](jtwpa-quantum-noise.md)
are available about the coherent pump or depleted multitone operating point.
Strong drives and commensurate-frequency aliases need
additional care and independent transient checks. The separate KIT/TLS and
resonant JPA noise APIs retain their documented scopes.
