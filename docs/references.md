# References and their role in the models

The papers below provide physical background and modeling context. A citation
does not mean that the code reproduces every published equation, geometry or
gain curve. The repository's finite-circuit harmonic balance, numerical
convergence and JoSIM comparisons are described in its own guides/reports.
Use the matching entries in [references.bib](../references.bib).

## KIT and dielectric loss

<a id="eom2012"></a>
**[eom2012]** B. H. Eom, P. K. Day, H. G. LeDuc and J. Zmuidzinas,
“A wideband, low-noise superconducting amplifier with high dynamic range,”
*Nature Physics* **8**, 623–627 (2012).
[DOI: 10.1038/nphys2356](https://doi.org/10.1038/nphys2356).
Background for nonlinear kinetic-inductance traveling-wave amplification.
The current KIT gain APIs here implement 3WM, rather than this paper's full
experimental 4WM device.

<a id="malnou2021"></a>
**[malnou2021]** M. Malnou, M. R. Vissers, J. D. Wheeler, J. Aumentado,
J. Hubmayr, J. N. Ullom and J. Gao, “Three-Wave Mixing Kinetic Inductance
Traveling-Wave Amplifier with Near-Quantum-Limited Noise Performance,”
*PRX Quantum* **2**, 010302 (2021).
[DOI: 10.1103/PRXQuantum.2.010302](https://doi.org/10.1103/PRXQuantum.2.010302).
Context for DC-biased KIT 3WM and signal/idler noise accounting. Example
preset values in this repository are not asserted to be that device.

<a id="martinis2005"></a>
**[martinis2005]** J. M. Martinis et al., “Decoherence in Josephson Qubits
from Dielectric Loss,” *Physical Review Letters* **95**, 210503 (2005).
[DOI: 10.1103/PhysRevLett.95.210503](https://doi.org/10.1103/PhysRevLett.95.210503).
Motivates resonant TLS dielectric loss and its drive dependence. It does not
calibrate the particular loss tangent or pump-hole linewidth used here.

<a id="gao2008"></a>
**[gao2008]** J. Gao et al., “Experimental evidence for a surface distribution
of two-level systems in superconducting lithographed microwave resonators,”
*Applied Physics Letters* **92**, 152505 (2008).
[DOI: 10.1063/1.2906373](https://doi.org/10.1063/1.2906373).
Background for surface participation and electric-energy weighting. The IDC
half-gap boundary problem and capacitance-ratio wiring remain separately
documented implementation choices.

## Resonant JPA

<a id="eichler2014"></a>
**[eichler2014]** C. Eichler and A. Wallraff, “Controlling the dynamic range
of a Josephson parametric amplifier,” *EPJ Quantum Technology* **1**, 2 (2014).
[DOI: 10.1140/epjqt2](https://doi.org/10.1140/epjqt2).
Principal reference for the driven Kerr-resonator/input-output description,
gain, noise and dynamic-range context. See [the JPA guide](resonant-jpa.md).

## Josephson traveling-wave circuits

<a id="obrien2014"></a>
**[obrien2014]** K. O'Brien, C. Macklin, I. Siddiqi and X. Zhang,
“Resonant Phase Matching of Josephson Junction Traveling Wave Parametric
Amplifiers,” *Physical Review Letters* **113**, 157001 (2014).
[DOI: 10.1103/PhysRevLett.113.157001](https://doi.org/10.1103/PhysRevLett.113.157001).
Reference for capacitively coupled LC phase matching in 4WM JJ lines.
`obrien_2014()` supplies published circuit values; the finite-port solver
does not claim to reproduce every continuum gain curve in the paper.

<a id="zorin2016"></a>
**[zorin2016]** A. B. Zorin, “Josephson Traveling-Wave Parametric Amplifier
with Three-Wave Mixing,” *Physical Review Applied* **6**, 034006 (2016).
[DOI: 10.1103/PhysRevApplied.6.034006](https://doi.org/10.1103/PhysRevApplied.6.034006).
Physical context for flux-biased one-junction rf-SQUID chains and tunable
quadratic/Kerr nonlinearities. The implementation uses its documented
unique-branch domain and finite-circuit current law.

<a id="frattini2017"></a>
**[frattini2017]** N. E. Frattini, U. Vool, S. Shankar, A. Narla,
K. M. Sliwa and M. H. Devoret, “3-wave mixing Josephson dipole element,”
*Applied Physics Letters* **110**, 222603 (2017).
[DOI: 10.1063/1.4984142](https://doi.org/10.1063/1.4984142).
Reference for the SNAIL element and its asymmetric potential. A small bare
quartic term is not treated here as proof of universally Kerr-free behavior.

<a id="roudsari2023"></a>
**[roudsari2023]** A. Fadavi Roudsari et al., “Three-wave mixing traveling-wave
parametric amplifier with periodic variation of the circuit parameters,”
*Applied Physics Letters* **122**, 052601 (2023).
[DOI: 10.1063/5.0127690](https://doi.org/10.1063/5.0127690).
Experimental SNAIL-TWPA and dispersion-engineering context. The synthetic
identically loaded examples here do not reproduce its periodic pattern.

## Quantum limits and distributed noise

<a id="caves1982"></a>
**[caves1982]** C. M. Caves, “Quantum limits on noise in linear amplifiers,”
*Physical Review D* **26**, 1817–1839 (1982).
[DOI: 10.1103/PhysRevD.26.1817](https://doi.org/10.1103/PhysRevD.26.1817).
Foundation for phase-preserving added-noise limits and the distinction from
phase-sensitive amplification. The code uses the finite-gain bound and
vacuum variance 1/2, not a universal half-quantum floor at every gain.

<a id="caves2012"></a>
**[caves2012]** C. M. Caves, J. Combes, Z. Jiang and S. Pandey,
“Quantum limits on phase-preserving linear amplifiers,”
*Physical Review A* **86**, 063802 (2012).
[DOI: 10.1103/PhysRevA.86.063802](https://doi.org/10.1103/PhysRevA.86.063802).
Broader phase-preserving noise-channel context. The present code uses
Gaussian second moments, rather than arbitrary noise distributions.

<a id="peng2022"></a>
**[peng2022]** K. Peng, M. Naghiloo, J. Wang, G. D. Cunningham, Y. Ye and
K. P. O'Brien, “Floquet-Mode Traveling-Wave Parametric Amplifiers,”
*PRX Quantum* **3**, 020306 (2022).
[DOI: 10.1103/PRXQuantum.3.020306](https://doi.org/10.1103/PRXQuantum.3.020306).
Context for Floquet channels, loss ports, backward waves and quantum
efficiency. The finite-circuit implementation here is separately derived
from its circuit Jacobian and resistor-bath normalization.

## Circuit simulator

<a id="josim2019"></a>
**[josim2019]** J. A. Delport, K. Jackman, P. le Roux and C. J. Fourie,
“JoSIM—Superconductor SPICE Simulator,” *IEEE Transactions on Applied
Superconductivity* **29**(5), 1300905 (2019).
[DOI: 10.1109/TASC.2019.2897312](https://doi.org/10.1109/TASC.2019.2897312).
JoSIM provides the independent classical RCSJ transient calculations.
See its [official documentation](https://joeydelp.github.io/JoSIM/) and
the repository's [backend guide](josim-backend.md). Optional fork-specific
correctors and per-run versions are recorded in the validation reports;
the simulator paper is not an attestation of those builds or quantum noise.

## Cite the software and the relevant physics

<a id="openems"></a>
**[openems]** Thorsten Liebig and openEMS contributors, *openEMS*.
[Official documentation](https://docs.openems.de/),
[ports and post-processing](https://docs.openems.de/en/latest/concepts/ports.html),
[CSXCAD polygon primitives](https://docs.openems.de/en/latest/python/CSXCAD/CSPrimitives/CSPrimPolygon.html).
The optional EM worker uses CSXCAD geometry and openEMS EC-FDTD port
admittance for capacitor extraction. These external projects have their own
licenses. Their documentation does not qualify this repository's illustrative
material stack or replace numerical convergence checks.

<a id="klayout"></a>
**[klayout]** Matthias Köfferlein and KLayout contributors, *KLayout*.
[Official project](https://www.klayout.de/),
[standalone Python API](https://www.klayout.org/klayout-pypi/) and
[DRC runsets](https://www.klayout.de/doc/manual/drc_runsets.html).
The optional layout package uses KLayout for polygon construction,
connectivity geometry, GDS exchange and DRC. KLayout is a separate dependency
with its own license; its documentation does not qualify the example processes.

Use [CITATION.cff](../CITATION.cff) for the software and record the exact
version/commit and numerical controls used. Add the references appropriate
to your device and noise model from [references.bib](../references.bib).
Synthetic benchmark outputs should be identified as calculations, not
independent experimental measurements.

## Inductance extraction

<a id="fasthenry1994"></a>
**[fasthenry1994]** M. Kamon, M. J. Tsuk and J. K. White,
“FASTHENRY: A Multipole-Accelerated 3-D Inductance Extraction Program,”
*IEEE Transactions on Microwave Theory and Techniques* **42**(9),
1750–1758 (1994). [DOI](https://doi.org/10.1109/22.310584).
The [author-hosted earlier DAC 1993 version](https://www.rle.mit.edu/cpg/publications/pub53.pdf)
describes volume filaments, partial-inductance matrices and terminal-network
reduction. The superconducting extension and this repository's finite-mask
adapter require separate verification; the current bounded benchmarks use
direct field evaluation, not multipole acceleration.

<a id="whiteley_fasthenry"></a>
**[whiteley_fasthenry]** S. R. Whiteley, *Xic Inductance/Resistance Extraction
Interface*: [geometry](https://wrcad.com/manual/xicmanual/node498.html),
[terminal definition](https://wrcad.com/manual/xicmanual/node499.html), and
[superconducting output extensions](https://wrcad.com/manual/xicmanual/node501.html).
Primary documentation for finite electrode equivalencing, ordered terminal
pairs and the Whiteley `lambda` extension. The upstream normal-metal solver
does not supply that extension. See the [independent London implementation](magnetic-extraction.md)
and its [native analytic and circuit checks](magnetic-jpa-validation.md).
