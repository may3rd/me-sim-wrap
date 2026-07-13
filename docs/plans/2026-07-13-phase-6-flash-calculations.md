# Phase 6 Flash Calculations Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add deterministic Rachford-Rice, PR TP, bubble-pressure, dew-pressure, and PH calculations for the frozen five-compound domain.

**Architecture:** Keep flash algorithms in one standard-library module using the Phase 5 `PengRobinsonMixture` implementation. Validate inputs before iteration, use bounded solves, and return structured convergence records instead of leaking partial states. Calculate total phase enthalpy from the existing ideal correlations plus PR departure enthalpy using one explicit reference state.

**Tech Stack:** Python 3 standard library, immutable dataclasses, `unittest`, existing compound/interaction/correlation JSON, DWSIM 9.0.4 golden cases.

---

## Equation references

- Rachford, H. H. and Rice, J. D., “Procedure for Use of Electronic Digital Computers in Calculating Flash Vaporization Hydrocarbon Equilibrium,” 1952.
- Michelsen, M. L., “The Isothermal Flash Problem. Part I. Stability,” *Fluid Phase Equilibria* 9, 1982.
- Michelsen, M. L., “The Isothermal Flash Problem. Part II. Phase-Split Calculation,” *Fluid Phase Equilibria* 9, 1982.

The implementation tests must record the exact edition, page/equation, or DOI used for each independent vector.

## Preconditions

Do not start implementation until the Phase 5 DWSIM parity capture passes and `codex/phase-5-pr` is merged. The flash solver consumes PR fugacity and departure properties directly.

Accepted Phase 6 scope:

- Rachford-Rice phase split
- TP flash
- isothermal bubble pressure
- isothermal dew pressure
- PH flash

PS is deferred until the Phase 10 expander requires it. TV and PV remain deferred until a unit operation requires them.

## Task 1: Solver result contracts and Rachford-Rice

**Files:**

- Create: `src/mesim/thermo/flash.py`
- Create: `tests/test_flash.py`

1. Write failing tests for:
   - all-liquid classification where `F(0) < 0`
   - all-vapor classification where `F(1) > 0`
   - a two-phase split with a known vapor fraction
   - `K = 1` degeneracy
   - invalid composition, nonpositive `K`, nonfinite input, and mismatched array lengths
2. Run `python3 -m unittest tests.test_flash -v` and confirm the expected failures.
3. Add immutable `SolverReport` and `FlashResult` dataclasses. Every result must include `converged`, `iterations`, `residual`, `algorithm`, `warnings`, `failure_reason`, and the last calculated state.
4. Implement the standard Rachford-Rice residual:

   ```text
   F(beta) = sum(z_i (K_i - 1) / (1 + beta (K_i - 1)))
   ```

5. Use bisection on `[0, 1]`; do not add SciPy. Stop when both the interval and residual tolerances pass or the iteration limit is reached.
6. Calculate phase compositions from the accepted vapor fraction and verify each composition is nonnegative and sums to one within `1e-12` absolute.
7. Return a structured failure on iteration exhaustion. Raise `ValidationError` only for invalid caller input.
8. Run the focused test, `python3 scripts/validate.py --quiet`, and the full suite.
9. Commit `feat: add bounded Rachford-Rice solver`.

## Task 2: PR phase-stability test

**Files:**

- Modify: `src/mesim/thermo/flash.py`
- Modify: `tests/test_flash.py`

1. Write failing vapor-like and liquid-like tangent-plane-distance tests using independent published or DWSIM-backed vectors.
2. Add Wilson initial estimates only for initialization:

   ```text
   ln(K_i) = ln(Pc_i / P) + 5.373 (1 + omega_i) (1 - Tc_i / T)
   ```

3. Implement both vapor-like and liquid-like trial phases. Iterate trial composition using PR fugacity coefficients and calculate tangent-plane distance.
4. Treat a negative minimum tangent-plane distance below tolerance as unstable. Do not use cubic-root count as a mixture stability test.
5. Bound trial iterations and return the last trial composition, residual, iteration count, and failure reason.
6. Test stable liquid, stable vapor, unstable feed, near-critical feed, zero mole fractions, and deterministic repeatability.
7. Run focused and full verification.
8. Commit `feat: add PR phase-stability test`.

## Task 3: TP flash

**Files:**

- Modify: `src/mesim/thermo/flash.py`
- Modify: `tests/test_flash.py`

1. Write failing TP cases for stable liquid, stable vapor, two-phase, near-critical, and zero-fraction components.
2. Run the focused test and confirm failure before implementation.
3. Validate temperature, pressure, ordered compound IDs, composition, and required interaction data.
4. Run the stability test first. Return one phase only when stability passes; do not infer single phase from Wilson `K` values or cubic-root count.
5. For an unstable feed:
   - initialize `K` with Wilson
   - solve Rachford-Rice
   - calculate liquid and vapor PR fugacity coefficients
   - update `ln(K_i) = ln(phi_i_liquid) - ln(phi_i_vapor)`
   - repeat until fugacity and material-balance residuals pass
6. Use log-space `K` values. If residuals oscillate or increase, halve the update damping down to `1/16`. If the bounded iteration still fails, return a structured failure; do not accept the last state as converged.
7. Require all of these convergence checks:
   - maximum component fugacity residual `<= 1e-8`
   - compound material-balance residual `<= 1e-10`
   - phase compositions sum to one within `1e-12` absolute
   - vapor fraction remains in `[0, 1]`
8. Run focused tests, quiet validation, and the full suite.
9. Commit `feat: add PR TP flash`.

## Task 4: Bubble and dew pressure

**Files:**

- Modify: `src/mesim/thermo/flash.py`
- Modify: `tests/test_flash.py`

1. Implement isothermal pressure calculations only; temperature bubble/dew calculations are deferred.
2. Write failing pure-component and mixture bubble/dew pressure tests.
3. Use a positive caller-supplied pressure bracket. Do not silently invent or widen engineering bounds.
4. At each pressure, iterate PR liquid/vapor fugacity coefficients before evaluating:
   - bubble residual `sum(z_i K_i) - 1`
   - dew residual `sum(z_i / K_i) - 1`
5. Solve pressure by bisection. Raise `ValidationError` for an invalid caller-supplied bracket; return a structured failure for iteration exhaustion.
6. Verify bubble pressure is not below dew pressure for the same valid two-phase envelope input without an explicit explanation.
7. Run focused and full verification.
8. Commit `feat: add PR bubble and dew pressure`.

## Task 5: Total phase enthalpy and PH flash

**Files:**

- Modify: `src/mesim/thermo/flash.py`
- Modify: `tests/test_flash.py`

1. Fix the calculation reference explicitly at `298.15 K` and `101325 Pa`; ideal enthalpy is zero at that reference for each compound.
2. Write failing tests for ideal-mixture enthalpy, phase enthalpy, total two-phase enthalpy, and a PH round trip generated from a converged TP state.
3. Calculate phase molar enthalpy as:

   ```text
   h_phase = sum(x_i * delta_h_ideal_i(T, 298.15 K)) + h_departure_PR
   ```

4. Reject temperatures outside any required heat-capacity correlation range. Do not extrapolate during flash calculations.
5. Require the caller to supply a positive finite temperature bracket for PH. Evaluate a complete TP flash at each trial temperature.
6. Solve temperature by bisection. The result converges only when both the enthalpy residual and the inner TP flash converge.
7. Include outer and inner iteration counts, final enthalpy residual in `J/kmol`, temperature bracket, warnings, and failure reason in the result.
8. Test single-liquid, single-vapor, phase-crossing, unreachable target enthalpy, invalid bracket, and deterministic repeatability.
9. Run focused and full verification.
10. Commit `feat: add bounded PR PH flash`.

## Task 6: DWSIM parity and exit gate

**Files:**

- Modify: `scripts/capture_dwsim_reference.ps1` only if existing generic property capture is insufficient
- Create: `tests/golden/pr-flash.json` on Windows
- Create: `tests/golden/pr-flash-repeat.json` on Windows
- Modify: `tests/test_flash.py`
- Modify: `docs/compatibility.md`

1. Capture DWSIM cases for stable vapor, stable liquid, two-phase TP, near-critical TP, bubble pressure, dew pressure, single-phase PH, and phase-crossing PH.
2. Capture twice and require `python scripts/validate.py --compare ...` to report normalized equality with no property read errors.
3. Compare temperature, pressure, and vapor fraction to DWSIM within `1e-6` relative unless a reviewed compound-data or reference-state difference is documented.
4. Independently verify compound material balance, phase-composition sums, fugacity equality, and PH energy closure; DWSIM agreement does not replace these invariants.
5. Run:

   ```bash
   python3 -m unittest tests.test_flash -v
   python3 -m unittest discover -s tests -v
   python3 scripts/validate.py --quiet
   git diff --check
   ```

6. Mark Phase 6 supported in `docs/compatibility.md` only after every gate passes.
7. Commit `test: verify PR flash calculations against DWSIM`.

## Exit gate

Phase 6 passes only when:

- every supported calculation converges or returns a structured failure
- no unconverged state is exposed as usable output
- TP component material balances close within `1e-10` relative
- converged two-phase states satisfy component fugacity equality within `1e-8`
- PH results close enthalpy within `max(1e-6 * abs(target), 1e-3) J/kmol`
- repeated local calculations are byte-identical
- captured DWSIM results meet the stated tolerance or carry a reviewed model-difference record
