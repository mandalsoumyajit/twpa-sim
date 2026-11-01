# Bug Tracker

Systematic audit and resolution log for defects identified in `twpa-sim`.

## Summary Table

| ID | Component | Location | Severity | Status | Verification |
| --- | --- | --- | --- | --- | --- |
| **B1** | EM / openEMS mesh generation | [`src/twpa/em/openems_worker.py`](../src/twpa/em/openems_worker.py) | High | `FIXED` | `tests.test_em` (`test_worker_mesh_rounding_avoids_float_drift_collision`) |
| **B2** | Kerr/Flux JPA steady-state solver | [`src/twpa/jpa.py`](../src/twpa/jpa.py) | High | `FIXED` | `tests.test_jpa` (`test_fold_root_reported_as_marginal`) |
| **B3** | Layout validation, `from_circuit` snapping & KIT turns | [`src/twpa/layout/bindings.py`](../src/twpa/layout/bindings.py), [`src/twpa/layout/generate.py`](../src/twpa/layout/generate.py), [`src/twpa/layout/kit.py`](../src/twpa/layout/kit.py) | High | `FIXED` | `tests.test_layout` (`test_electrical_sizing_and_incompatible_capacitance`) |
| **B4** | IDC Laplace solver mesh resolution & electrostatic validation | [`src/twpa/idc.py`](../src/twpa/idc.py), [`src/twpa/electrostatics.py`](../src/twpa/electrostatics.py) | High | `FIXED` | `tests.test_solvers` (`test_electrostatics_parallel_plate_and_dielectrics`) |
| **B5** | Physical SQUID JPA flux validation & mode root bracket | [`src/twpa/jpa_circuit.py`](../src/twpa/jpa_circuit.py) | Medium-High | `FIXED` | `tests.test_josim` (`test_squid_equilibrium_and_bias_ramp`) |
| **B6** | Coupled JPA / LayoutNetwork degeneracy & Port 0 reference | [`src/twpa/em/coupled_jpa.py`](../src/twpa/em/coupled_jpa.py), [`src/twpa/em/network.py`](../src/twpa/em/network.py) | Medium | `FIXED` | `tests.test_magnetic`, `tests.test_layout_annotation` |
| **B7** | KIT nonlinear `route_profile`, `compression` & `nonlinear` API | [`src/twpa/nonlinear.py`](../src/twpa/nonlinear.py) | Medium | `FIXED` | `tests.test_pipeline` (`test_nonlinear_small_signal_and_gradient_zero`) |
| **B8** | 4WM JTWPA cubic nonlinearity analytical Jacobian | [`src/twpa/jtwpa.py`](../src/twpa/jtwpa.py) | Low-Medium | `FIXED` | `tests.test_jtwpa` (`test_sparse_jacobian`) |
| **B9** | Validation, provenance & layout extraction edge cases | [`src/twpa/models.py`](../src/twpa/models.py), [`src/twpa/em/magnetic_annotation.py`](../src/twpa/em/magnetic_annotation.py), [`src/twpa/em/backannotation.py`](../src/twpa/em/backannotation.py), [`src/twpa/em/fasthenry.py`](../src/twpa/em/fasthenry.py), [`src/twpa/em/magnetic.py`](../src/twpa/em/magnetic.py), [`src/twpa/layout/validate.py`](../src/twpa/layout/validate.py) | Low-Medium | `FIXED` | `tests.test_solvers`, `tests.test_em`, `tests.test_magnetic` |

---

## Detailed Bug Reports

### B1: Floating-point mesh line collision in `openems_worker.py`
- **Status**: `FIXED`
- **Location**: [`src/twpa/em/openems_worker.py`](../src/twpa/em/openems_worker.py)
- **Root Cause**: `np.arange(low[axis], high[axis]+dx/2, dx)` accumulates IEEE 754 floating-point drift (`~1.78e-14`) for non-dyadic decimal step sizes (`mesh_um = 0.1, 0.15, 0.2, 0.3`). When combined with exact polygon vertices, `-half`, `half`, and `0`, `np.unique` fails to merge algebraically identical coordinates, triggering `np.min(np.diff(lines)) < dx*1e-5` (`ValueError: Near-duplicate mesh lines; choose a compatible mesh step`).
- **Fix**: Snap generated mesh lines to floating-point precision (`np.round(..., 12)`) before `np.unique` so algebraically coincident coordinates merge cleanly while genuine sub-threshold geometric collisions still raise `ValueError`.
- **Verification**: Verified via `tests.test_em.EMContract.test_worker_mesh_rounding_avoids_float_drift_collision` (8/8 tests passing).

### B2: Premature fold-root tolerance in `pump_points` and `select_pump_point` branch jump
- **Status**: `FIXED`
- **Location**: [`src/twpa/jpa.py`](../src/twpa/jpa.py)
- **Root Cause**: `pump_points` evaluates `abs(polynomial(x)) <= tol` at stationary turning points `turns` (`tol = 64 * eps * drive`) and skips adjacent brackets `[lo, hi]` whenever either endpoint has `abs(f) <= tol`. Because the cubic response is quadratic near a fold ($\Delta P \propto (\Delta x)^2$), drives strictly below the upper fold by up to `tol` have two distinct roots separated by $\Delta x \sim \sqrt{\text{tol}} \sim 10^{-7}$ (one stable lower-branch root and one unstable middle-branch root), both of which are skipped in favor of a single `'marginal'` root at `turn0`. `select_pump_point(..., branch='lower')` then ignores the lower branch and jumps to the upper branch. Drives slightly above the fold also record a spurious `'marginal'` root at `turn0`.
- **Fix**: Compute `photon_flux = power/(h*frequency)` directly to avoid `abs(sqrt(...) * exp(1j*phase))**2` rounding drift, tighten fold knot tolerance to `4*np.finfo(float).eps*max(drive, np.finfo(float).tiny)`, and tighten root deduplication to `1e-11*max(1., abs(root))`.
- **Verification**: Verified via `tests.test_jpa.JPATests.test_fold_root_reported_as_marginal` (14/14 tests passing).

### B3: Layout validation vs. generation off-grid half/quarter dimensions and KIT turn bracket validation
- **Status**: `FIXED`
- **Location**: [`src/twpa/layout/bindings.py`](../src/twpa/layout/bindings.py), [`src/twpa/layout/generate.py`](../src/twpa/layout/generate.py), [`src/twpa/layout/kit.py`](../src/twpa/layout/kit.py)
- **Root Cause**:
  1. `JosephsonLayout.validate` and `KITLayout.validate` only check `p.coordinate(v)` on full dimensions, whereas `_josephson`, `Canvas.jj`, `Canvas.capacitor`, and `_kit` place geometry at half-dimensions (`v/2`) and quarter-dimensions (`capacitor_height * 0.5 / 2` for traveling-wave end half-shunt capacitors).
  2. `from_circuit` in `bindings.py` snaps `junction_width` and `capacitor_height` to `process.grid` instead of `2 * process.grid` (or `4 * process.grid` for 4WM `Chain`), producing validated layouts that fail in `generate()`.
  3. `_kit` computes `pw = min(d.pad_width, length/3)` and passes `-pw/2` to `c.box` without snapping `pw` to `2 * c.p.grid`.
  4. `KITLayout.validate` does not verify that the turn spine width equation in `kit.py` is bracketed on `[w_line_straight/10, 10*w_line_straight]` or that inner turn fingers do not cross $r \le 0$.
- **Fix**: Validate half/quarter dimensions in `JosephsonLayout.validate` and `KITLayout.validate`, snap `from_circuit` dimensions to the required grid multiples, snap `pw` in `_kit` to `2 * c.p.grid`, and validate KIT turn geometry/bracketing in `KITLayout.validate` and `kit.py`.
- **Verification**: Verified via `tests.test_layout.LayoutTests.test_electrical_sizing_and_incompatible_capacitance` (11/11 tests passing).

### B4: Unresolved thin oxide / electrode in `idc_fields` and `electrostatic` oxide mask validation
- **Status**: `FIXED`
- **Location**: [`src/twpa/idc.py`](../src/twpa/idc.py), [`src/twpa/electrostatics.py`](../src/twpa/electrostatics.py)
- **Root Cause**:
  1. In `idc_fields`, when `0 < oxide_thickness < mesh/2`, no cell centers `(xc, yc)` fall inside `[0, oxide_thickness]`, silently returning `participation = 0.0`.
  2. In `idc_fields`, when `mesh` is too coarse to place grid nodes on the metal thickness `[0, metal_thickness]` or across the gap (`ix >= nx - 1`), `metal` is empty or fills the domain, yielding `v = 0` and `participation = nan`.
  3. `electrostatic` accepts non-boolean `oxide` arrays without validation.
- **Fix**: Validate in `idc_fields` that the mesh resolves the gap, metal electrode (`len(jm) >= 1`), and conformal oxide layer when `oxide_thickness > 0` (`np.any(oxide)`), raising a descriptive `ValueError` when the mesh is too coarse. Validate in `electrostatic` that `oxide` is boolean-valued (or 0/1) when provided.
- **Verification**: Verified via `tests.test_solvers.Solvers.test_electrostatics_parallel_plate_and_dielectrics` (8/8 tests passing).

### B5: `SQUIDJPA.junction_phase` flux validation and `SQUIDJPA._mode` root-bracketing range
- **Status**: `FIXED`
- **Location**: [`src/twpa/jpa_circuit.py`](../src/twpa/jpa_circuit.py)
- **Root Cause**:
  1. `SQUIDJPA.junction_phase(self, flux=None)` does not validate `flux` when explicitly supplied (`abs(flux) < 0.45` and finite).
  2. `SQUIDJPA._mode` searches only `[guess*0.8, guess*1.2]`, failing when `arm_inductance` and `junction_capacitance` shift the fundamental resonance outside $\pm 20\%$ of the lumped LC guess or when `guess*1.2` crosses the series-arm pole $\omega_{\text{pole}} = \sqrt{(L_J + L_{\text{arm}})/(L_J L_{\text{arm}} C_J)}$.
- **Fix**: Validate explicit `flux` arguments in `junction_phase` and widen the fundamental resonance bracket in `_mode` while staying strictly below the junction/arm pole.
- **Verification**: Verified via `tests.test_josim.CircuitMeasurementTests.test_squid_equilibrium_and_bias_ramp` (11/11 tests passing).

### B6: Missing degeneracy checks in `coupled_jpa.response` and `LayoutNetwork.response`, and `LayoutNetwork` Port 0 validation
- **Status**: `FIXED`
- **Location**: [`src/twpa/em/coupled_jpa.py`](../src/twpa/em/coupled_jpa.py), [`src/twpa/em/network.py`](../src/twpa/em/network.py)
- **Root Cause**:
  1. `coupled_jpa.response` does not reject `f == point.frequency` (degenerate 4WM where $f_s = f_i = f_p$).
  2. `LayoutNetwork.response` only checks `np.any(abs(fs) < 1.)`, missing 3WM degeneracy at $f = 0.5 f_p$ and $f = 1.5 f_p$ where `abs(fs)` contains colliding frequencies (`+0.5*fp` and `-0.5*fp`).
  3. `LayoutNetwork.__post_init__` does not require non-empty `ports` and `junctions`, and `LayoutNetwork.netlist` references `SRC0` to `'0'` while `port_current(0)` uses `p['negative']`.
- **Fix**: Reject degenerate/colliding sideband frequencies in both `coupled_jpa.response` and `LayoutNetwork.response`, validate non-empty `ports` and `junctions` in `LayoutNetwork.__post_init__`, and connect `SRC0` to `node(self.data['ports'][0]['negative'])` in `LayoutNetwork.netlist`.
- **Verification**: Verified via `tests.test_magnetic` and `tests.test_layout_annotation` (all 21 tests passing).

### B7: `route_profile('ideal')` return type, `compression` crossing logic, and `nonlinear` options/scalar validation
- **Status**: `FIXED`
- **Location**: [`src/twpa/nonlinear.py`](../src/twpa/nonlinear.py)
- **Root Cause**:
  1. `route_profile('ideal')` returns scalar `0.` instead of `0.*s`, breaking shape consistency when evaluated on a NumPy array `s`.
  2. `compression` checks `if len(crossing) and crossing[0] > 0:`, which fails to detect a valid upward crossing `(deficit[:-1] < 1) & (deficit[1:] >= 1)` if `deficit[0] >= 1`.
  3. `nonlinear` does not accept `**pump_options` to forward to `select_pump(ref, bias, **pump_options)` and does not validate that `signal_frequency` is a scalar float.
- **Fix**: Return `0.*s` in `route_profile('ideal')`, detect the first upward crossing `np.flatnonzero((deficit[:-1] < 1) & (deficit[1:] >= 1))` in `compression`, validate scalar `signal_frequency` in `nonlinear`, and forward `**pump_options` to `select_pump`.
- **Verification**: Verified via `tests.test_pipeline` (6/6 tests passing).

### B8: Finite-difference Jacobian inaccuracy in `jtwpa._nonlinear(..., nonlinearity='cubic', derivative=True)`
- **Status**: `FIXED`
- **Location**: [`src/twpa/jtwpa.py`](../src/twpa/jtwpa.py)
- **Root Cause**: `_nonlinear` uses a centered finite difference (`step = 1e-6`) when `nonlinearity == 'cubic'` and `derivative=True`, whereas `nonlinearity == 'sine'` uses the exact analytical Jacobian. The finite difference introduces $O(\text{step}^2) \sim 10^{-13}$ truncation error even at `branch = 0`.
- **Fix**: Replace the finite-difference approximation in `_nonlinear(..., nonlinearity='cubic', derivative=True)` with the exact analytical Wirtinger derivatives of `_cube(branch) / 6`.
- **Verification**: Verified via `tests.test_jtwpa.JTWPA.test_sparse_jacobian` (9/9 tests passing).

### B9: Input validation, provenance hashing, and layout inductor terminal sorting edge cases
- **Status**: `FIXED`
- **Location**: [`src/twpa/models.py`](../src/twpa/models.py), [`src/twpa/em/magnetic_annotation.py`](../src/twpa/em/magnetic_annotation.py), [`src/twpa/em/backannotation.py`](../src/twpa/em/backannotation.py), [`src/twpa/em/fasthenry.py`](../src/twpa/em/fasthenry.py), [`src/twpa/em/magnetic.py`](../src/twpa/em/magnetic.py), [`src/twpa/layout/validate.py`](../src/twpa/layout/validate.py)
- **Root Cause**:
  1. `positive()` in `models.py` accepts `True` because `isinstance(True, int)` is `True` in Python.
  2. `magnetic_annotation.certify_magnetic` allows `primary` to be passed inside `controls={'frequency': primary}` when `primary` is at `>= 1e11` Hz, producing zero frequency sensitivity (`df = 0`).
  3. `backannotation.back_annotate_jpa` silently ignores `inductive_coupler=True` in a JPA manifest instead of directing callers to `back_annotate_network_jpa`.
  4. `fasthenry.probe` passes `self.executable` directly to `sha256sum` / `Path(self.executable)`, which fails if `executable` is a command name on `PATH` rather than a file path.
  5. `magnetic.layout_magnetic_job` contains a duplicate `if patch_side_um > min(b[2]-b[0], b[3]-b[1]):` check.
  6. `layout/validate.py` orders `IND` marker cut terminals by the global bounding-box center `pieces[k].bbox().center().x` of the entire `M1` island rather than the local contact with the inductor cut `cr`.
- **Fix**: Reject booleans in `positive()`, `Section`, and `Design`, require `fcontrol` frequency to differ from `primary` frequency in `magnetic_annotation.py`, reject `inductive_coupler=True` in `back_annotate_jpa`, resolve `self.executable` via `shutil.which` / `which` in `fasthenry.py`, remove the duplicate check in `magnetic.py`, and sort `IND` terminals in `layout/validate.py` by their local contact center `(pieces[k] & cr.sized(1)).bbox().center().x`.
- **Verification**: Verified via `tests.test_solvers`, `tests.test_em`, and `tests.test_magnetic` (all passing).
