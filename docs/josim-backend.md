# JoSIM backend

The optional `twpa.backends.JoSIM` runner executes physical Josephson-junction
netlists and reads CSV traces. It supports a native JoSIM executable or a Linux
executable reached from Windows through WSL. JoSIM is installed separately;
no binary or simulator source is bundled with this package.

## Configure and run

```python
from twpa.backends import JoSIM
from twpa.jpa_circuit import ParallelJPA, Tone, Transient, jpa_netlist

simulator = JoSIM(executable="josim-cli")
# On Windows, an explicit WSL executable can be used instead:
# simulator = JoSIM(executable="/path/to/josim-cli", wsl_distribution="Ubuntu")
circuit = ParallelJPA(critical_current=13e-6, shunt_capacitance=999e-12)
deck = jpa_netlist(circuit, [Tone(frequency=1e9, power=1e-22)],
                   Transient(step=.5e-12, stop=3e-6, start=1e-6))
result = simulator.run(deck, "validation/runs/example")
print(result.traces.keys())
```

The run directory must be new. It contains the exact input deck, output CSV,
stdout/stderr logs and a manifest with simulator version, command arguments,
netlist SHA256, metadata, timing and return code. Failures and timeouts preserve
these artifacts. Commands are argument lists without shell interpolation.
Windows paths with spaces are supported. WSL runs also require GNU `timeout`
inside Linux so a deadline kills the simulator rather than only its launcher.

The runner uses JoSIM's default **phase formulation** (`analysis=1`). Voltage
formulation is available with `analysis=0`. In the tested JoSIM build,
voltage-mode coupled-loop runs accumulated a small static flux-offset error
that became larger when the timestep was halved; this affected quadrature axes.
The phase formulation is checked in the validation report. Do not infer that a
smaller timestep alone guarantees correctness in every simulator formulation.

`corrector=N` requests an optional predictor-corrector feature of the enhanced
JoSIM build. The runner probes `-h` first and rejects this request if unsupported.
The default is zero, using the simulator's ordinary algorithm. CSV is used for
upstream compatibility rather than relying on fork-specific binary formats.

## Physical resonator circuits

`ParallelJPA` is an RCSJ junction in parallel with a shunt capacitor and internal
resistance, driven through a Thevenin source resistor. The source resistance
is also the port reference impedance. Junction capacitance appears **once** in
the JJ model; `shunt_capacitance` is additional. `RTYPE=0` and equal RN/R0 give
explicit constant junction damping, without gap switching or automatic thermal
noise. Temperature-dependent model scaling is not requested.

Its reduced model uses the classical LC frequency, total capacitance C, external
decay `1/(Z0*C)`, and internal conductance divided by C. The weak-phase Kerr
coefficient is `-e^2/(2*C*hbar)`. This classical circuit mapping does not add a
vacuum normal-ordering frequency correction or calibrate a fabricated device.

`SQUIDJPA` contains two identical RCSJ junctions, explicit loop-arm inductors,
a driven coil and equal/opposite mutual couplings. The magnetic drive is a coil
current, not a fictitious JJ flux terminal. The bias current ramps from zero:
starting a voltage-mode simulation with a pre-existing DC coil current did not
establish the intended loop flux in the tested build. Keep the ramp and discard
subsequent settling time even when using phase analysis.

The stable zero-fluxoid branch solves

    chi + (L_arm*Ic/phi_bar)*sin(chi) = pi*Phi_external/Phi0

Linearized arm/junction admittances determine the unpumped resonance, effective
capacitance (half the susceptance slope), and loss. The effective flux pair
drive comes from the derivative of that resonance with applied flux. It is not
fitted to the pumped gain. Only the supported stable bias range and weak flux
modulation are intended; screening, parasitic modes and higher drive require
additional convergence/domain checks.

For sinusoidal coil flux `Phi_ac*sin(omega_p*t+phase)`, the voltage-envelope
convention gives `lambda = -0.5*(d omega_mode/d flux)*Phi_ac*exp(-i*phase)`.
Flux amplitudes here are normalized by Phi0. Frequencies are Hz, lambda and
decay rates are angular rates. The effective flux JPA frame is half the physical
pump frequency. A single-mode rotating-wave approximation is used only by the
reduced model; JoSIM integrates the RCSJ circuit itself.

## Measure signals at the port

`twpa.measurement.port_waves` uses instantaneous or peak-phasor voltage/current
with current directed into the circuit:

    V_inc = (V + Z0*I)/2
    V_ref = (V - Z0*I)/2
    power = abs(V_wave)^2/(2*Z0)

Source voltage is twice the incident-wave amplitude. Comparing tank voltage
with source voltage directly would give a different, incorrect gain definition.
`fit_tones` jointly fits DC, cosines and sines using actual CSV timestamps and
returns peak phasors under the `exp(-i*omega*t)` convention. It checks sampling,
tone separation and fit conditioning.

`sideband_measurement` subtracts a pump-only run at identical timestamps and
fits signal, idler and pump harmonics together. Idler conversion is returned in
photon-flux units, including the signal/idler frequency ratio. Multiple probes
must be sufficiently weak to make their mutual intermodulation negligible.
`degenerate_measurement` instead uses two orthogonal probe phases to identify
the real two-quadrature transfer matrix. It does not pretend that coincident
signal and idler are independent modes.

## Reproduce the validation

```powershell
python tools/validate_jpa_josim.py --executable /path/to/josim-cli `
  --wsl Ubuntu --output validation/runs/jpa-0p5ps --step 5e-13 --compression
```

Omit `--wsl` for native execution. Add `--corrector 4` only for a simulator that
advertises this feature. The script records model and circuit parameters,
measurements, errors and explicit pass/fail criteria. Output folders are never
overwritten. Compare 1 ps and 0.5 ps steps and repeat a longer settling/window
run; the carrier period alone is not an adequate convergence criterion for a
narrow resonance. These values are for the supplied approximately 1 GHz test
circuits, not universal simulation settings.

An optional live integration test is enabled by `TWPA_JOSIM_EXECUTABLE`, with
`TWPA_JOSIM_WSL` for Windows/WSL and `TWPA_JOSIM_TEST_ROOT` to select a writable
temporary directory. Ordinary CI uses analytic/backend tests and stored circuit
regressions and does not require JoSIM to be installed.

The [validation results](josim-validation.md) distinguish classical quantities
verified by JoSIM from quantum predictions. JoSIM's classical resistor noise
does not validate vacuum fluctuations, output squeezing noise or added photons.
The physical circuit tests validate gain and deterministic quadrature response;
quantum noise remains tested through input-output identities and analytic limits.

Simulator references: [JoSIM documentation](https://joeydelp.github.io/JoSIM/)
and [thermal-noise implementation](https://joeydelp.github.io/JoSIM/tech_disc/).

## Four-wave traveling-wave circuits

`twpa.jtwpa.chain_netlist` generates finite RCSJ chains with physical source
and load resistors, half shunts at the ends, and optional capacitively coupled
LC resonators. `twpa.josim_jtwpa.run_chain` stores the shared circuit metadata;
`measure_chain` extracts two-port waves and optionally subtracts a pump-only
trajectory. `measure_bloch` extracts the interior voltage recurrence from
three neighboring nodes. `Tone.ramp_time` adds exponential source turn-on to
reduce long-lived startup ringing. See the [4WM model guide](four-wave-jtwpa.md)
and [JTWPA circuit comparisons](jtwpa-josim-validation.md) for API conventions,
actual benchmark results and reproduction commands.
