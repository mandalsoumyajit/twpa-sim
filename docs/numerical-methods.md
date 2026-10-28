# Understanding the numerical calculations

The examples use several numerical methods because they answer different
physical questions. This guide explains what the methods do and what their
diagnostics mean. It complements the [physical overview](physical-overview.md)
and the equations in the individual device guides.

## Passive circuits before pumping

A passive calculation asks how a sufficiently weak wave propagates when the
amplifier pump is absent. Linear circuit theory then applies. For a repeated
line cell, an **ABCD matrix** relates voltage and current at one end of the
cell to voltage and current at the other. Cascading these matrices gives the
finite line's reflection and transmission.

**Bloch analysis** instead asks how a wave changes from one identical cell to
the next in an infinite periodic line. The complex cell-to-cell multiplier
determines phase advance and attenuation. The finite circuit can still have
reflections and resonances even when its periodic cell is in a passband.
This is why the tutorials inspect both propagation and physical port response.

## Solving a pumped operating point

The strong pump changes the state about which a weak signal propagates.
In an effective JPA, a steady pump can have several possible cavity photon
occupations. The solver finds those roots and checks their linear stability.
Selecting a branch chooses a physical operating point, not just a numerical
root ordering. In this model, small perturbations must decay for the branch
to be stable.

The finite Josephson-chain solvers use **harmonic balance**. They represent
each node's periodic or multitone motion as a sum of retained sinusoids,
insert that representation into the nonlinear circuit equations, and solve
for the unknown complex amplitudes. Kirchhoff's current law, abbreviated
**KCL**, requires current entering each node to balance current leaving it.
The **residual** is the remaining mismatch in those equations after solving,
in the normalization specified by the solver.

A **Jacobian** is the matrix of derivatives of those equations with respect
to the unknown amplitudes. Newton's method uses it to improve a trial solution.
The weak-signal calculation also uses the Jacobian to determine how a small
signal perturbs the pumped circuit. **Continuation** means using a nearby
solution, for example at a slightly lower pump power, as the next initial
guess. It helps follow a solution branch through a parameter sweep.

The KIT propagation solvers make a different approximation. They evolve
slowly varying pump, signal and idler envelopes along the line rather than
solving every reflected wave at every circuit node. Their ordinary
differential equations, abbreviated **ODEs**, describe how those envelopes
change with propagation distance or cell index. The forward-only boundary
is part of that model, not merely a faster implementation of the finite
Josephson circuit.

## Harmonics, sidebands and integration grids

A nonlinear current creates frequencies that were absent from the input.
Retaining only pump, signal and idler may omit a harmonic or conversion path
that changes the predicted gain. `harmonics` and `sidebands` determine which
frequency components a finite-circuit solve keeps. The device guides give
their precise indexing; the same numerical value does not select identical
frequencies in the 3WM and 4WM APIs.

**Quadrature** in a numerical setting means a weighted sampling rule for an
integral. Harmonic balance samples the nonlinear current over drive phases
and projects it onto the retained Fourier components. Here `quadrature`
controls the number of phase samples. This use of the word differs from a
field's cosine/sine quadratures in a squeezing calculation.

Refine three things separately: the number of retained harmonics, the number
of sidebands, and the number of integration samples. A small residual only
shows that the equations for the chosen truncation were solved accurately.
It does not show that omitted frequency components are negligible.

## Convergence and stability answer different questions

**Convergence** asks whether a reported observable changes acceptably when
the numerical approximation is refined. The observable might be signal gain,
idler conversion, capacitance or mutual inductance. The acceptable change
must be chosen for that use case and compared with material uncertainty.

**Dynamical stability** asks whether a small perturbation grows in time.
A converged nonlinear root can be unstable. The effective resonant JPA tests
its drift-matrix eigenvalues. The extracted five-state JPA tests a periodic
pump using a **Floquet multiplier**: an eigenvalue describing perturbation
growth over one pump period. A multiplier magnitude below one means decay.
The distributed Josephson-chain harmonic-balance APIs do not currently
perform a general dynamical stability test.

The phase guard is another distinct check. It limits the calculated junction
or element phase excursion to the model's supported regime. Passing the guard
does not establish stability or predict the device's switching threshold.

## What the noise diagnostics check

The Gaussian noise calculation propagates small fluctuations about the chosen
classical state. A periodically driven circuit couples a fluctuation at one
frequency to other frequencies separated by pump-related intervals. These
coupled components are often called **Floquet channels**. Retaining more
channels can change the noise even when the central signal gain changes little;
[Peng et al.](references.md#peng2022) discuss the physical role of these modes.

An output field must retain the canonical quantum commutation relations.
The reported **commutator error** measures the defect in that algebra after
including both external inputs and internal loss baths. The covariance must
also satisfy the uncertainty relation. These are strong consistency checks,
but neither forces harmonic convergence or validates an assumed loss model.
The [noise guide](jtwpa-quantum-noise.md) defines the matrices and the
[Caves references](references.md#caves1982) give the physical noise limit.

## Comparing with a transient circuit simulation

JoSIM evolves a classical physical circuit in time. To measure weak-signal
gain in the presence of a much larger pump, the validation scripts compare a
pump-plus-probe run with a pump-only run and fit the difference at the signal
and idler frequencies. This separates coherent conversion from the pump's
large background and generated harmonics.

The time step must resolve relevant circuit oscillations. The discarded
startup time must allow transients to settle. The retained acquisition window
must distinguish nearby fitted tones: its duration needs to be long compared
with the inverse of their frequency separation. Refine these controls
independently. JoSIM validates classical circuit response; it does not supply
an independent quantum-vacuum noise calculation.

## Refining a layout extraction

A field mesh divides a conductor or dielectric into numerical elements.
Refine the mesh without changing the physical geometry or terminal footprints.
For electric extraction, also enlarge the surrounding domain and extend the
simulation duration. For London inductance, refine the in-plane current basis,
the film-thickness subdivision and the basis near contacts. The
[magnetic workflow](magnetic-extraction.md) describes its additional kinetic
energy and signed-coupling checks.

Convergence estimates numerical sensitivity within the selected approximation.
It does not determine an unknown penetration depth, qualify fabrication rules,
or show that an omitted package or electromagnetic mode is unimportant.
