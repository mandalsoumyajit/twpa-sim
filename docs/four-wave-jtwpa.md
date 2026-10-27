# Four-wave Josephson traveling-wave amplifiers

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

c = jtwpa.Chain(critical_current=PHIBAR/2.5e-9,
                ground_capacitance=1e-12, junction_capacitance=1e-12,
                cells=384)
c = replace(c, port_impedance=float(jtwpa.bloch(c, [1e9]).impedance[0].real))
power = .5*c.port_impedance*(.4*c.critical_current)**2
r = jtwpa.four_wave_response(c, 1e9, power, .9e9, power*1e-6)
print(10*np.log10(r.gain), r.idler_conversion, r.pump_transmission)
```

Run `python examples/four_wave_jtwpa.py` to compare uniform and RPM chains.
`linear_response` returns exact finite S11/S21 using scaled ABCD cascades.
`bloch` returns principal-zone wavenumber in rad/cell, attenuation in
nepers/cell, Bloch impedance and the corresponding lossless-cell passband
classification. Divide wavenumber by `cell_length` to obtain rad/m. Near poles
and band edges use finite scattering and adequately resolved frequency grids.

## Equations and normalization

For node flux Phi and branch phase delta=(Phi_left-Phi_right)/phi_bar,
phi_bar=hbar/(2e), the branch current is

    I = Ic*sin(delta) + Cj*phi_bar*delta_ddot + phi_bar*delta_dot/Rj.

Ground capacitors, optional ground resistors, and optional resonator
admittances enter node KCL. The model solves this finite circuit directly in
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

By default `nonlinearity='sine'` projects the exact sine current onto a
two-phase Fourier basis. Retained frequencies are

    f(h,m) = h*fp + m*(fs-fp),

with positive odd h through `harmonics` (default 3), and integer m from
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
