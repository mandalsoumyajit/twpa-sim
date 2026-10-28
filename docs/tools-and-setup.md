# Tools and setup

TWPA Sim supplies Python amplifier models, layout generators, extraction
adapters and comparison scripts. Some stages use separate engineering tools.
“Optional” below means optional for the package as a whole: a tool becomes
required when you choose the stage that uses it. Installing the core package
does not install all circuit or electromagnetic simulators.

If you want to explore gain and quantum noise from circuit parameters, start
with the Python environment alone. If you want a checked, layout-derived
resonant JPA model, follow the [design cycle](design-workflow.md) and add the
layout, field-extraction and circuit-comparison tools as you reach those stages.

## Which tools do I need?

| Tool | Role in the design cycle | Required when | How it is obtained |
| --- | --- | --- | --- |
| Python 3.10+, NumPy and SciPy | Solve electrical models, propagation, gain and quantum fluctuations; orchestrate other tools | Every Python simulation | Install Python, then `pip install -e .`; pip installs the declared NumPy/SciPy dependencies |
| Matplotlib | Make plots in scripts that offer plotting | You request those plots | `pip install -e '.[plot]'` |
| KLayout **Python module** | Build/read GDS, inspect polygons, check annotated connectivity and run the Python geometry/DRC operations | Layout generation/checking and preparing geometry for extraction | `pip install -e '.[layout]'` |
| KLayout **application/executable** | View masks and markers; independently run the repository's native Ruby design-rule deck | Visual inspection or that independent native DRC step | Separate [KLayout installation](https://www.klayout.de/build.html) |
| JoSIM executable | Integrate physical Josephson circuits in time and compare measured classical response with the Python model | Actual JoSIM circuit validation, rather than netlist generation/dry runs | Separate [JoSIM installation](https://joeydelp.github.io/JoSIM/#initial-setup) |
| openEMS **and CSXCAD Python bindings**, with their native libraries | Build the supported capacitor field model and extract capacitance | The openEMS electric-extraction stage | Separate [openEMS installation and Python-interface setup](https://docs.openems.de/en/latest/install/index.html) |
| Whiteley superconducting FastHenry executable | Extract signed self/mutual inductance including linear London kinetic energy | The JPA magnetic-interconnect extraction stage | Separate [Whiteley distribution/build](https://www.wrcad.com/freestuff.html), or a compatible selected installation |
| Linux `timeout`, `prlimit`, `env`, `sha256sum`, `which`; WSL on Windows for this runner | Bound and identify Linux jobs; reach Linux executables from Windows | The current native-Linux/WSL FastHenry adapter; WSL JoSIM also uses GNU `timeout` and comparison scripts resolve the selected binary | Linux tools in the selected system/distribution; WSL is installed separately on Windows |

The standalone KLayout Python package and full application are different
installations. The Python package includes geometry/report functionality and
does not provide the full application interface or native Ruby runner; see
the [official distinction](https://www.klayout.org/klayout-pypi/). Installing
`.[layout]` therefore enables the Python checks, while `--klayout` selects
the additional executable-based check.

CSXCAD supplies the structure, material and mesh description used by openEMS.
Both compiled Python interfaces must work in the interpreter selected for the
worker. AppCSXCAD, if installed with openEMS, can help inspect its structures;
the repository's numerical worker does not require that viewer. MATLAB and
Octave are not required by this Python workflow.

## Set up the main Python environment

From a source checkout on Windows:

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -e .
.\.venv\Scripts\python examples/howto/03_kerr_jpa.py
```

Use `.venv/bin/python` on Linux/macOS. Below, `python` means this main
interpreter, whether you activate the environment or call its full path.
Install additional Python features as needed:

```powershell
python -m pip install -e '.[layout,plot]'
python -c "import numpy, scipy, twpa; print('Electrical model imports work')"
python -c "import klayout.db; print('KLayout Python', klayout.db.__version__)"
```

These commands do not install JoSIM, the KLayout application, openEMS native
libraries or FastHenry. Their executables/interpreters are passed explicitly
to the relevant adapter. Git is useful for obtaining the checkout and recording
the revision; a C++ compiler and build tools are needed only if you choose
to build an external simulator rather than use a compatible binary.

## Check KLayout generation and independent DRC

With the Python layout extra installed, first generate/check the small
examples without a native application:

```powershell
python examples/layout/generate_all.py --output validation/runs/layout-python-check
```

For the additional native Ruby DRC step, install the full application, check
its executable and pass it explicitly:

```powershell
klayout -v
python examples/layout/generate_all.py --output validation/runs/layout-native-check --klayout klayout
```

On Windows the executable may be `klayout_app.exe`; use its quoted full path
for `--klayout` if needed. Inspect the JSON/marker reports and the native
`.lyrdb` files rather than treating a zero process return code as a clean
DRC report. The [layout guide](layout.md#independent-geometry-and-drc-decks)
explains which checks each path runs. A full application's version and its
standalone Python package version may differ, so record both.

## Add JoSIM for an independent circuit calculation

JoSIM reads a SPICE-like netlist of junctions, capacitors, inductors, resistors
and drives, then writes classical voltage/phase traces. TWPA Sim constructs
the netlist and measures incident/reflected/transmitted waves from those traces.
JoSIM does not extract fields from GDS and its classical noise output is not
a quantum-vacuum reference. See the [backend guide](josim-backend.md).

Install or build JoSIM following its [official instructions](https://joeydelp.github.io/JoSIM/).
Check the actual executable before launching a long study:

```powershell
josim-cli -v
josim-cli -h
```

Use the full executable path if it is not on `PATH`. From Windows, a Linux
installation can instead be checked with:

```powershell
wsl -d Ubuntu --exec /path/to/josim-cli -v
wsl -d Ubuntu --exec /path/to/josim-cli -h
```

Configure the Python runner with the same location:

```python
from twpa.backends import JoSIM

runner = JoSIM(executable='josim-cli')  # Native executable.
# Alternative for a Linux executable launched from Windows:
# runner = JoSIM('/path/to/josim-cli', wsl_distribution='Ubuntu')
print(runner.probe())
```

The ordinary backend defaults to CSV output and no enhanced corrector.
`--corrector 4` in benchmark commands is a capability of the selected enhanced
build, not a requirement satisfied by every JoSIM installation. Request it
only if the executable advertises the feature. Similarly, the London-network
validation driver defaults to native named-float64 `.npy` traces for long,
finely sampled records: that particular path requires a build advertising
NPY output. The runner checks these capabilities instead of assuming them.
Its `--trace-format csv` option still requires adequate timestamp precision;
it is not a guarantee that coarse-precision CSV can resolve that benchmark.

Keep the same executable throughout a convergence study and record its
version and binary hash. The repository's [validation reports](README.md#validation-and-provenance)
identify the builds actually used, rather than treating an arbitrary newer
version as the same numerical reference.

## Add openEMS and CSXCAD for electric extraction

Follow the [openEMS installation guide](https://docs.openems.de/en/latest/install/index.html),
including its Python-interface instructions and compiled-library requirements.
Use an interpreter compatible with those bindings. It can be separate from
the main environment running TWPA Sim; the worker needs NumPy, openEMS and
CSXCAD but does not need to import `twpa` or KLayout.

```text
Main Python: twpa + NumPy/SciPy + optional KLayout/Matplotlib
    writes geometry/material job.json
    starts the explicitly selected worker interpreter
openEMS Python: NumPy + compiled openEMS/CSXCAD interfaces and native libraries
    returns extraction.json, field-model inputs and logs
Main Python: checks convergence and back-annotates the circuit
```

As an installation check, run this code **in the openEMS interpreter**:

```python
import openEMS
import CSXCAD
print('openEMS', openEMS.__version__)
print('CSXCAD', CSXCAD.__version__)
```

On Windows, the native libraries may need their DLL directory registered
before either import. For a direct check, retain the directory handle:

```python
import os
library_directory = os.add_dll_directory(r'C:\path\to\openEMS')
import openEMS
import CSXCAD
print(openEMS.__version__, CSXCAD.__version__)
```

For an extraction run, pass the corresponding paths instead:
`--openems-python C:/path/to/openems/python.exe` and, on Windows when needed,
`--openems-dll-directory C:/path/to/openEMS`. The worker handles library
registration. The [electric guide](em-extraction.md) describes the actual MIM
electrodes, fit band, finite domain and independent refinement checks.
The current adapter extracts a capacitor component, not all chip/package
capacitances or nonlinear amplifier gain.

The repository also provides a runtime probe that uses the same import and
DLL-registration path as the extraction worker. In PowerShell:

```powershell
& C:/path/to/openems/python.exe src/twpa/em/openems_worker.py --probe --dll-directory C:/path/to/openEMS
```

Omit `--dll-directory` where it is unnecessary. This prints installed openEMS
and CSXCAD versions without launching a field simulation. An import failure
here is an environment/binary-library problem to resolve before changing
capacitor geometry or mesh settings.

## Add superconducting FastHenry for magnetic extraction

This flow needs the Whiteley superconducting extension, which supports London
penetration-depth specifications. Ordinary normal-metal MIT FastHenry input
does not implement that extension; the distinction is documented by
[Whiteley](https://wrcad.com/manual/xicmanual/node501.html). Select a compatible
installed executable and verify its banner and hash. The current adapter runs
it natively on Linux or through WSL from Windows; native Windows/macOS
FastHenry execution is not implemented by this adapter.

For a Windows/WSL setup, create a local installation file such as the following.
Replace all example paths and the hash with values for your installation:

```json
{
  "variant": "Whiteley superconducting FastHenry",
  "linux_binary": "/path/to/fasthenry",
  "wsl_distro": "Ubuntu",
  "binary_sha256": "REPLACE_WITH_THE_EXECUTABLE_SHA256"
}
```

For native Linux, use the same variant, `linux_binary` and `binary_sha256`;
`wsl_distro` is not used. Check the binary with `fasthenry -h` and
`sha256sum /path/to/fasthenry` on Linux, or prefix those commands with
`wsl -d Ubuntu --exec` from Windows. Then test the adapter without solving:

```python
from twpa.em import FastHenry
runner = FastHenry.from_installation('fasthenry-installation.json')
print(runner.probe())
```

Pass this file explicitly with `--installation`. This avoids depending on
another project's default installation location. The runner compares the
binary hash, requires the expected Whiteley banner, and uses GNU wall-time
and address-space limits. GNU coreutils `stdbuf` flushes diagnostic output
before a timeout; `timeout`, `prlimit` and `sha256sum` are also required.
Each default case is capped at 180 s and 4 GiB
address space. The runner explicitly bounds dense-library threads to two
by default (`FastHenry(..., threads=...)` accepts 1–4); the local accelerated
Whiteley build also honors its recorded dense-backend thread setting. The CLI supports one or two such
workers; these limits bound jobs, not the RAM needed by every possible design.

The [magnetic guide](magnetic-extraction.md) explains finite terminal patches,
physical returns and signed mutual coupling. Total London inductance already
contains kinetic energy; the separately reported energy decomposition is not
another inductance to add to the circuit.

## Confirm readiness for your chosen path

For an electrical study, the core imports and HOWTO examples are sufficient
to start. For layout, also run a small generation/check case. For independent
native DRC, verify the KLayout executable and inspect its actual report.
For field extraction and JoSIM, first check each selected runtime above,
then run the supported [JPA design cycle](design-workflow.md#worked-jpa-cycle).

External executables have their own installation requirements and licenses;
they are not redistributed inside TWPA Sim. Use explicit paths, record tool
versions, and preserve the layout/material snapshots and run inputs so the
next design iteration can be compared with the previous one.

The general [layout-network driver](layout-backannotation.md) has separate,
explicit 8 GiB/900 s starting limits and up to two workers. Its Cartesian
film mesh has a 40,000-filament guard. These bounds can reject larger designs;
the adapter does not silently replace a failed full network with a uniform
inductance estimate. Its total-L certificate does not include the dedicated
JPA's separate kinetic-energy decomposition. The KIT calibration example
requires KLayout and the Python core, with no native EM executable.
