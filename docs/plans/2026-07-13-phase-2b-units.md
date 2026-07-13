# Phase 2B Process Units Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend the explicit unit boundary with verified process units required through the internal-alpha scope.

**Architecture:** Keep the existing immutable unit table and conversion functions. Add exact aliases by semantic dimension; retain affine conversion only for absolute temperature and gauge pressure. No parser or new dependency.

**Tech Stack:** Python 3 standard library and `unittest`.

---

### Task 1: Correct molar-property bases

**Files:**
- Modify: `tests/test_units.py`
- Modify: `src/mesim/units.py`

1. Add failing vectors proving `1 kJ/kmol = 1000 J/kmol` and `1 J/mol = 1000 J/kmol`.
2. Run `python3 -m unittest tests.test_units -v` and confirm failure.
3. Make `J/kmol` the molar enthalpy SI base and correct its aliases.
4. Run focused and full tests.
5. Commit `fix: correct molar property unit bases`.

### Task 2: Add alpha process dimensions

**Files:**
- Modify: `tests/test_units.py`
- Modify: `src/mesim/units.py`

1. Add failing equation and round-trip vectors for temperature difference, gauge pressure, time, mass, amount, volumetric flow, mass/molar heat capacity, mass/molar entropy, specific volume, kinematic viscosity, diffusivity, surface tension, pressure gradient, mass flux, heat flux, acceleration, and force.
2. Run the focused test and confirm unknown-unit failures.
3. Add the minimum explicit aliases required by those vectors.
4. Run focused tests, `python3 scripts/validate.py --quiet`, full tests, and `git diff --check`.
5. Commit `feat: add explicit process unit dimensions`.

### Task 3: Correct the Phase 2 gate

**Files:**
- Modify: `docs/plans/2026-07-13-python-process-simulator-rewrite.md`
- Modify: `docs/compatibility.md`

1. Record Phase 2B coverage and state that full DWSIM display-unit parity is not required.
2. Make future unit additions demand-driven by verified model inputs and outputs.
3. Run the full suite.
4. Commit `docs: define process unit coverage boundary`.

