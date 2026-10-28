# Interconnect, kinetic inductance and coupling extraction

Single-junction resonant JPA layouts now support a hybrid field-to-circuit
flow: openEMS extracts the actual shunt-capacitor electrodes, and a
London-aware Whiteley FastHenry build extracts the complete finite metal
current network, including return paths, contacts and signed mutual coupling.
The annotated full-sine circuit is compared with JoSIM without fitting a
linewidth or gain. See the [numerical report](magnetic-jpa-validation.md).

## Materials and reference planes

Install `twpa-sim[layout]` to prepare jobs. The separate
[magnetic stack JSON](../src/twpa/em/stacks/london_illustrative.json) specifies
film thicknesses, dielectric gap, penetration depths and their temperature.
The example assumes 200 nm films and 90 nm penetration depth at 20 mK;
these are illustrative inputs, not a measurement or a SeeQC process model.
Supply measured operating-temperature penetration depths for a real design.
No temperature extrapolation, vortex or nonlinear kinetic model is implicit.
The [electric stack](../src/twpa/em/stacks/mim_illustrative.json) and magnetic
stack must have the same physical dielectric gap.

`jpa_magnetic_job` reads the actual M1, M2 and VIA polygons, checks the layout
against its manifest and design rules, and tiles the whole finite metal area.
The film current basis connects tile centers to shared tile faces. All mask,
contact and terminal boundaries are explicit cuts. Both in-plane directions
are present; conductor areas must equal the mask areas. Contacts include the
film-center-to-interface paths and the physical interlayer fill. Shared film
paths are not duplicated when a terminal coincides with a contact.

Every terminal is a fixed finite electrode. Its nodes are connected only
within that electrode using `.equiv`; there is no idealized ground plane or
shortcut across a return conductor. JJ and capacitor reference planes are
the facing film surfaces. Vertical half-film currents retain London energy.
The feed field inside each locally ideal electrode is excluded by this
partition. Mesh refinement keeps the electrode footprint fixed. Comparing
another footprint changes the reference-plane model; the present convergence
certificate does **not** establish footprint independence.

The generator supports two actual wiring topologies:

| Layout | Raw oriented conductor ports | Circuit currents |
| --- | --- | --- |
| `inductive_coupler=true` | CAP_M2→JJ_M1, JJ_M2→CAP_M1, coil_lo→coil_hi | Closed JJ/capacitor tank and finite pickup loop |
| Original bare JPA | EXT→JJ_M1, EXT→CAP_M2, JJ_M2→GROUND, CAP_M1→GROUND | JJ feed and capacitor feed, both including the common physical return |

For the pickup case, `I_raw = T I_loop` with
`T = [[1,0],[1,0],[0,1]]`; for the bare case
`T = [[1,0],[0,1],[1,0],[0,1]]`. The extraction retains the native signed
complex port matrix and reduces it as `Z_loop = T.T @ Z_raw @ T`.
Only then is `L = Im(Z_loop)/(2*pi*f)` formed. Reversing a port reverses its
mutual entries. A negative mutual inductance is valid. Nonreciprocal,
nonfinite, incomplete or nonpositive-energy matrices are rejected.

## London energy and avoiding double counting

The Whiteley build uses `lambda` in the geometry length unit (um here).
With positive penetration depth and `sigma=0`, the model is lossless London
magnetoquasistatics. Its inductance matrix already includes magnetic and
kinetic energy. **Use that total matrix once in the circuit.** Do not add an
extra sheet-kinetic inductance afterward.

For fixed terminal currents, the discrete lossless problem minimizes

`U(s) = U_magnetic[J] + s U_kinetic[J]`,

where all penetration depths are multiplied by `sqrt(s)`. At the stationary
solution, the first-order change from redistribution of `J` vanishes, so

`L_kinetic = dL/ds at s=1`, and `L_magnetic = L_total - L_kinetic`.

The implementation evaluates centered differences at `s=1±0.1` and
`s=1±0.05`, compares them, and checks nonnegative kinetic energy. This is an
energy decomposition for the selected fixed current basis and linear London
model, not a separate fitted series element. A uniform-wire native control
checks `Lk = mu0*lambda^2*length/(width*thickness)` and lambda-squared scaling.

`L_SC - L_PEC` is reported separately: it also includes changes in magnetic
energy caused by current redistribution and is **not** the pure kinetic
matrix. Exact `lambda=0,sigma=0` is not a valid native conductor. The PEC
comparison is a positive-small-lambda limit, checked at factors 1e-3 and 1e-4.

## Run the flow

The [commented example](../examples/em/london_jpa.py) prepares an audited
pickup-loop GDS and input deck without starting a native solver:

```powershell
python examples/em/london_jpa.py --output validation/runs/london-example
```

For a full comparison, install the **Whiteley superconducting** FastHenry
variant, openEMS/CSXCAD in a separate Python environment, and JoSIM. Standard
normal-metal FastHenry does not supply the required London model. External
simulators retain their own licenses; no binaries or vendor PDK tables are
redistributed here. Runtime and geometry methods follow the finite-electrode,
explicit-return flow reviewed in `spec2si-seeqc`; this is an independent
adapter using public illustrative material inputs.

```powershell
python tools/validate_magnetic_jpa.py --output validation/runs/london-full `
  --installation C:/path/fasthenry-installation.json `
  --openems-python C:/path/openems/python.exe `
  --openems-dll-directory C:/path/openems `
  --executable /path/to/josim-cli --wsl Ubuntu --corrector 4 `
  --step 2.5e-15
```

Omit `--wsl` for native JoSIM. Use `--corrector` only with an executable that
advertises it. The default trace format is precision-preserving named float64
NPY; the selected JoSIM build must advertise `.npy` output. `--trace-format csv`
uses CSV at a coarser output interval but still rejects duplicate timestamps.
No timestamp repair is performed. Most other JoSIM APIs retain CSV by default.

The installation JSON has `variant: "Whiteley superconducting FastHenry"`,
`linux_binary`, `binary_sha256`, and, on Windows, `wsl_distro`. Instead of a
configuration file, instantiate `FastHenry(executable, wsl_distribution,
binary_sha256)` directly. Native Linux runs use GNU `timeout` and `prlimit`.
The runner pins the binary SHA, uses single-thread BLAS, and imposes 180 s
and 4 GiB address-space limits per case. The CLI permits one or two concurrent
cases; use `--workers 1` on smaller machines. Dense direct LU is selected for
a single frequency; iterative direct-field solves serve frequency sweeps.
The local mesh guard rejects jobs above 20,000 specified filaments.

The driver creates fresh directories, saves decks, job/material/source hashes,
native matrices and solver logs, and requires these independent controls:

- Bulk pitch, film-thickness discretization and local contact basis: below
  2% change in each self inductance and 5% in signed mutual inductance.
- Independent 100 GHz lossless-London frequency control at fixed current basis.
- Converged PEC limit and kinetic-energy derivative, without changing geometry.
- Electric mesh, domain and duration controls from the [openEMS flow](em-extraction.md).
- Passive complex S11, pumped signal and idler against JoSIM, and extra pump
  harmonics/Floquet sidebands. Refine timestep and settling/acquisition as in
  the [validation report](magnetic-jpa-validation.md).

`--gds` accepts an existing supported single-JJ layout with its manifest.
`--extract-magnetic-only` stops after the magnetic certificate. To reuse
completed extraction, pass `--magnetic` and `--electric` certified artifacts;
they must match the **same GDS bytes, manifest, process and design**. Reusing a
similar-looking or newly regenerated file with another hash is rejected.

```python
from twpa.em import back_annotate_network_jpa, load_circuit

circuit, provenance = back_annotate_network_jpa(
    'jpa.gds', 'electric-certified.json', 'magnetic-certified.json',
    output='backannotation.json')
# CoupledJPA for a pickup loop; GalvanicJPA for the original signal/ground port.
assert circuit == load_circuit('backannotation.json')
```

The pickup network uses a 50 ohm external termination by default; the bare
galvanic example defaults to 5 kohm. These are explicit embedding assumptions,
not extracted resistances. Internal and junction loss resistances are explicit
too. `Cj` comes from junction area and the process density and is added once.
`Cs` comes from the actual capacitor electrodes; each attaches to its actual
internal circuit node. The full matrix becomes two inductors and a signed
JoSIM `K` element. Python uses the same five-state full-sine circuit, periodic
harmonic balance, Floquet stability and sideband response.

## Scope

This implements interconnect/return/contact inductance, linear kinetic energy
and inductive pickup or shared-return coupling for the supported resonant
JPA geometry. It does not yet extract electric stray/cross capacitance outside
the MIM coupon, packages, radiation, dielectric loss, bias coils, SQUID/SNAIL
loops or TWPA unit cells. The result is a finite-conductor MQS network plus a
component electric extraction. A flat London L versus frequency is a model
consistency check, not evidence that omitted full-wave effects are negligible
at 100 GHz. Recheck the model boundary and reference planes for another layout.

See [FastHenry and Whiteley references](references.md#fasthenry1994), especially
the [terminal definition](https://wrcad.com/manual/xicmanual/node499.html)
and [geometry construction](https://wrcad.com/manual/xicmanual/node498.html).
