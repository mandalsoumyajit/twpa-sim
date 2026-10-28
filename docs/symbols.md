# Symbols and technical terms

This page is a lookup aid for the device guides and examples. Subscripts are
written with underscores in code-style equations: `f_p` means pump frequency.
The guides define symbols locally as well. In particular, `lambda` is a
penetration **length** in the London extraction and a coupling **rate** in an
effective flux-pumped resonator; these are different physical quantities.

## Frequencies, waves and measurements

| Symbol or term | Meaning | Units or convention |
| --- | --- | --- |
| `f`, `f_p`, `f_s`, `f_i` | Ordinary frequency; pump, signal and idler frequencies | Hz |
| `omega`, `omega_0` | Angular frequency and unpumped angular resonance; `omega=2*pi*f` | rad/s |
| `t` | Time | s |
| `P`, `P_p` | Power and available incident pump power | W; available power specifies an incident wave, not absorbed power |
| `R`, `Z_0` | Real port reference resistance/impedance | ohm |
| `G`, `G_dB` | Signal power gain and `10*log10(G)` | Dimensionless; dB |
| `S11`, `S21` | Complex input reflection and forward transmission coefficients | Squared magnitudes give power ratios for compatible real reference ports |
| `k` | Bloch phase advance per repeated cell in the propagation guides | rad/cell; divide by cell length for rad/m |
| Passband / stopband | Frequency interval supporting / suppressing propagation in the corresponding lossless periodic cell | A finite lossy circuit still needs its scattering calculation |
| Harmonic | An integer multiple of a drive frequency | For example, the third pump harmonic is `3*f_p` |
| Sideband | A generated frequency offset from a carrier or pump harmonic | Its index is a numerical label, not a separate fabricated component |
| Phasor | Complex amplitude encoding the size and phase of a sinusoid | The guide states whether it multiplies `exp(+i*omega*t)` or `exp(-i*omega*t)` |

## Junctions, loops and resonators

| Symbol or API field | Meaning | Units or convention |
| --- | --- | --- |
| `I_c`, `critical_current` | Critical current in the stated physical junction | A; for SNAIL this refers to a large junction |
| `I_dc`, `I_scale` | Applied direct current and characteristic current for the KIT inductance law | A |
| `delta`, `delta_0`, `x` | Junction/element phase, equilibrium phase, and displacement from equilibrium | rad |
| `Phi`, `Phi_ext` | Node flux or applied magnetic flux, as specified locally | Wb; node-flux time derivative is voltage |
| `Phi_0`, `phi_bar` / `PHIBAR` | Flux quantum `h/(2*e)` and reduced flux quantum `Phi_0/(2*pi)` | Wb |
| `h`, `hbar`, `e` | Planck constant, `h/(2*pi)`, positive elementary charge | J s, J s, C |
| `L_J`, `L_loop` | Small-signal Josephson inductance and geometric loop inductance | H |
| `C_j`, `C_g`, `C_s` | Junction capacitance, node-to-ground capacitance, and shunt capacitance | F; local SNAIL formulas may use `C_s` for the small junction instead |
| `R_j` | Junction loss resistance in an RCSJ circuit | ohm |
| `E_J`, `U` | Josephson energy scale `phi_bar*I_c` and circuit potential energy | J |
| `beta`, `beta_L` | rf-SQUID screening parameter `L_loop*I_c/phi_bar` | Dimensionless |
| `alpha`, `n` in SNAIL | Small/large critical-current ratio and number of large series junctions | Dimensionless ratio; integer count |
| `c2`, `c3`, `c4` | Derivatives of `U/E_J` with respect to phase at equilibrium | Dimensionless; control linear curvature, 3WM and quartic nonlinearity |
| `a`, `a_in`, `a_out` in JPA equations | Intracavity amplitude, incoming and outgoing wave amplitudes | `abs(a)**2` is cavity photons; wave squares are photons/s |
| `kappa`, `kappa_external`, `kappa_internal` | Total energy-decay rate and its port/internal contributions | rad/s in the API; amplitude decays at `kappa/2` |
| `K`, `kerr` | Resonance angular-frequency shift per intracavity photon in the effective JPA model | rad/s per photon |
| `Delta` | Angular drive detuning from resonance; sign defined locally | rad/s |

## Loss and quantum noise

| Symbol or term | Meaning | Units or convention |
| --- | --- | --- |
| `T`, `k_B` | Bath temperature and Boltzmann constant | K; J/K |
| `n(f,T)` | Thermal mean photon occupation at frequency `f` | Dimensionless; zero for a zero-temperature bath |
| `nu=n+1/2` | Symmetrized thermal noise including vacuum | Photon units per spectral bin |
| `N_added` | Input-referred added noise | Quanta in the stated symmetrized convention |
| `a_dagger` | Creation operator, the Hermitian conjugate of annihilation operator `a` | Quantum notation; not an additional classical tone |
| `U`, `V`, `Ub`, `Vb` in noise matrices | External annihilation/creation transfer matrices and internal-bath counterparts | Photon-flux normalized; `U` here is not potential energy |
| `C`, `M` in noise matrices | Symmetrized ordinary and anomalous output covariances | `C` here is a covariance, not a circuit capacitance |
| `I`, `Omega` in matrix identities | Identity matrix and canonical quadrature commutator matrix | `I` here is not electrical current |
| TLS, participation | Two-level systems in a dielectric; fraction of electric energy stored in the lossy region | Participation is dimensionless |
| Loss tangent | Ratio characterizing dielectric dissipation relative to energy storage | Dimensionless |

## Layout and field extraction

| Symbol or term | Meaning | Units or convention |
| --- | --- | --- |
| M1, M2, VIA | First/second conducting layers and a physical interlayer contact | Layer roles in the process profile |
| GDS, DRC, manifest | Layout exchange file, design-rule checking, and the accompanying design/process/provenance record | DRC is a geometry check, not electrical performance validation |
| MIM capacitor | Metal–insulator–metal capacitor | Its two electrodes are separated by a dielectric |
| Reference plane / terminal patch | Physical boundary at which a field-extracted element connects to its circuit | Moving it changes which fields belong to the extracted element |
| `Z`, `Y` | Impedance and admittance matrices relating terminal currents and voltages | ohm; siemens |
| `L`, `M` in magnetic extraction | Self/total inductance matrix and signed mutual inductance | H; `M` here is not quantum anomalous covariance |
| `lambda`, `sigma` | London penetration depth and normal conductivity in the extraction model | um in the emitted geometry; conductivity follows the solver convention |
| `mu_0` | Vacuum permeability | H/m |
| `J`, `s` in the energy decomposition | Current density and an artificial multiplier of penetration-depth squared | A/m²; dimensionless |
| PEC, MQS, FDTD | Perfect electric conductor, magnetoquasistatic approximation, finite-difference time-domain method | Modeling assumptions/methods, not process qualifications |

For the physical picture, begin with the [overview](physical-overview.md).
For exact API normalization and mode ordering, use [conventions](conventions.md).
