# Python Process Simulator Rewrite Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a deterministic, cross-platform Python process-simulation kernel that grows from a verified hydrocarbon subset into functional parity for DWSIM thermodynamics, compounds, steady-state unit operations, reactions, columns, hydraulics, and dynamics.

**Architecture:** Use a Python calculation kernel with no dependency on the DWSIM runtime. Store compound and interaction data as immutable, versioned records; calculate internally in SI while preserving every submitted value and unit at the API boundary; expose the kernel through FastAPI only after direct Python calculations pass. Start with five light hydrocarbons, Peng-Robinson, TP/PH flash, and six basic unit operations. Add models only after the previous phase passes equation-level, DWSIM-reference, and conservation checks.

**Tech Stack:** Python 3.12, standard library first, NumPy and SciPy when numerical scope requires them, FastAPI and Pydantic at the HTTP boundary, `unittest` for runnable checks, JSON for versioned data and golden cases, Docker Linux/amd64 and Linux/arm64.

---

## 1. Non-negotiable boundaries

1. The new engine must run natively on macOS arm64 and Linux amd64/arm64.
2. The new package must never import, invoke, or load a DWSIM assembly.
3. DWSIM remains a Windows-only reference oracle until parity is accepted; it is not a production dependency.
4. Initial domain: methane, ethane, propane, n-butane, and n-pentane with Peng-Robinson vapor-liquid equilibrium.
5. Initial mode: deterministic steady state. Dynamics starts only after steady-state flowsheet convergence is stable.
6. Every input and output carries an explicit unit. Preserve the submitted value and unit; store a separate SI value for calculation. Never silently round or repair data.
7. Every compound constant, correlation, and binary interaction parameter carries source, source version, original unit, and import timestamp.
8. Solver failures are values, not crashes: return algorithm, iteration count, residuals, last state, and failure reason.
9. Python is the implementation language. C/C++ is permitted only for a profiled hotspot with a benchmark proving the need and a Python reference implementation proving equivalence.
10. “All DWSIM functionality” means the compatibility matrix in section 5. Desktop editors, icons, reports, spreadsheets, and CAPE-OPEN hosting are not calculation kernels and are outside the rewrite unless separately approved.

## 2. Repository target

```text
me-sim-wrap/
├── data/
│   ├── compounds/v1.json
│   ├── interactions/pr-v1.json
│   └── reactions/
├── docs/
│   ├── compatibility.md
│   └── plans/
├── scripts/
│   ├── capture_dwsim_reference.ps1
│   └── validate.py
├── src/mesim/
│   ├── __init__.py
│   ├── errors.py
│   ├── units.py
│   ├── compounds.py
│   ├── streams.py
│   ├── flowsheet.py
│   ├── thermo/
│   │   ├── ideal.py
│   │   ├── peng_robinson.py
│   │   ├── flash.py
│   │   ├── activity.py
│   │   ├── steam.py
│   │   └── electrolyte.py
│   ├── unitops/
│   │   ├── basic.py
│   │   ├── pressure.py
│   │   ├── heat.py
│   │   ├── separation.py
│   │   ├── hydraulics.py
│   │   ├── reactors.py
│   │   ├── columns.py
│   │   └── dynamics.py
│   └── api.py
├── tests/
│   ├── golden/
│   ├── test_units.py
│   ├── test_compounds.py
│   ├── test_thermo.py
│   ├── test_flash.py
│   ├── test_unitops.py
│   └── test_flowsheets.py
└── pyproject.toml
```

Create files only when their phase begins. Do not scaffold empty future modules.

## 3. Stable calculation contracts

These contracts are introduced only when the first implementation needs them. Do not create abstract base classes with one implementation.

### Quantity at trust boundaries

```python
@dataclass(frozen=True, slots=True)
class Quantity:
    value: float
    unit: str
    si_value: float
```

`value` and `unit` are immutable source data. `si_value` is derived explicitly by `convert_to_si()` and must never overwrite them.

### Compound record

```python
@dataclass(frozen=True, slots=True)
class Compound:
    id: str
    name: str
    cas: str
    formula: str
    molecular_weight_kg_per_kmol: float
    critical_temperature_k: float
    critical_pressure_pa: float
    acentric_factor: float
    normal_boiling_point_k: float
    provenance: dict[str, str]
```

Add optional properties only when a calculation consumes them. Missing required data raises `MissingCompoundData`; it is never estimated silently.

### Stream calculation state

```python
@dataclass(frozen=True, slots=True)
class StreamState:
    temperature_k: float
    pressure_pa: float
    molar_flow_kmol_s: float
    composition: tuple[float, ...]
    compound_ids: tuple[str, ...]
    enthalpy_j_per_kmol: float | None = None
    vapor_fraction: float | None = None
```

Validate finite values, positive absolute temperature and pressure, aligned compound arrays, nonnegative flow and composition, and composition sum within the configured solver tolerance.

### Result and failure

Every iterative result includes `converged`, `iterations`, `residual`, `algorithm`, and `warnings`. A non-converged result must not be applied to a stream or flowsheet.

## 4. Verification policy

Every phase must pass all applicable gates before the next phase starts.

### Gate A — equation checks

- Compare equations against a cited primary reference or standard.
- Check dimensional consistency.
- Test limiting cases, invalid inputs, and phase/root selection.
- Test deterministic repeatability using exact serialized inputs.

### Gate B — reference cases

- Generate golden JSON on x64 Windows using the vendored DWSIM revision.
- Record DWSIM commit/revision, property package, flash algorithm, compound data source, input values with units, outputs with units, and solver messages.
- Never update expected results because the Python result changed. A golden update requires a source/version note and review.

### Gate C — conservation and invariants

- Material balance per compound.
- Total mass and energy balance.
- Mole fractions remain within `[0, 1]` and sum to one within tolerance.
- Temperatures, pressures, densities, heat capacities, and transport properties stay inside model validity ranges.
- Repeated execution returns byte-identical JSON except explicit timestamps.

### Initial acceptance tolerances

| Check | Acceptance |
|---|---:|
| Unit conversion | Exact for scale-only conversions; `1e-12` relative otherwise |
| Composition sum | `1e-12` absolute |
| EOS compressibility and fugacity | `1e-8` relative to independent equation cases |
| DWSIM flash temperature, pressure, vapor fraction | `1e-6` relative or documented model difference |
| Compound material balance | `1e-10` relative |
| Unit-operation energy balance | `1e-8` relative |
| Full flowsheet outputs | `1e-5` relative after convergence |

Tolerances are versioned configuration used by validation scripts, not scattered constants.

### Required commands

```bash
python -m unittest discover -s tests -v
python scripts/validate.py
python scripts/validate.py --quiet
```

`validate.py` exits `0` only when every enabled compatibility case passes and `1` otherwise.

## 5. Compatibility matrix and order

### Thermodynamic models

| Phase | Scope |
|---|---|
| T0 | Ideal gas, ideal liquid, pure-component vapor pressure and heat-capacity correlations |
| T1 | Peng-Robinson and PR78, vapor and liquid roots, fugacity, departure enthalpy/entropy |
| T2 | TP, PH, bubble-pressure and dew-pressure flashes; stability testing. Add PS with the expander and TV/PV only on demand |
| T3 | Soave-Redlich-Kwong, Lee-Kesler-Plocker, Chao-Seader, Grayson-Streed |
| T4 | Wilson, NRTL, UNIQUAC, UNIFAC, modified UNIFAC; VLE, LLE, VLLE |
| T5 | IAPWS-IF97 steam/water and seawater |
| T6 | Electrolyte ideal, electrolyte SVLE, LIQUAC-style scope, sour water |
| T7 | Black oil, petroleum characterization, solids/SLE, hydrates, PC-SAFT and specialty EOS |
| T8 | Multi-phase Gibbs minimization and remaining flash variants |

Models are added because an approved flowsheet requires them, not merely because DWSIM contains them.

### Unit operations

| Phase | Scope |
|---|---|
| U0 | Material stream, energy stream, mixer, splitter, heater, cooler, valve, equilibrium separator |
| U1 | Pump, compressor, expander, component separator |
| U2 | Two-stream heat exchanger and shell-and-tube thermal rating |
| U3 | Pipe, fittings, orifice plate, relief valve; single- and two-phase hydraulic correlations |
| U4 | Conversion, equilibrium and Gibbs reactors |
| U5 | CSTR and PFR with kinetic reactions |
| U6 | Shortcut column, absorber and reboiled absorber |
| U7 | Rigorous staged columns and column solver variants |
| U8 | Recycle, energy recycle, adjust, specification and flowsheet convergence |
| U9 | Dynamic vessels, dynamic heat exchangers, controllers and time integration |
| U10 | Solar, wind, hydroelectric, fuel-cell and electrolyzer models |

### Explicit non-core parity items

- Spreadsheet unit operation: replace with typed import/export, not Excel automation.
- Python-script unit operation: add a restricted callback API only when a real internal model needs it.
- CAPE-OPEN host: separate adapter project; never couple it to the calculation kernel.
- UI editors, gauges, switches and reports: web-application work, not engine rewrite work.

## 6. Implementation phases

### Phase 0: Freeze reference and scope

**Files:**
- Create: `docs/compatibility.md`
- Create: `scripts/capture_dwsim_reference.ps1`
- Create: `tests/golden/schema.json`

**Steps:**

1. Record the vendored DWSIM revision and GPLv3 provenance in `docs/compatibility.md`.
2. Inventory property packages, flash algorithms, compound sources, unit operations, and DWSIM sample flowsheets.
3. Define the golden-case JSON schema with immutable input values and units.
4. Write the Windows capture script to load a case through `Automation3`, solve it, and emit normalized JSON without rounding.
5. Capture five compound records and these initial cases: single vapor, single liquid, two-phase flash, mixer, heater, valve, separator.
6. Run the capture twice and compare hashes; stop if outputs are nondeterministic.
7. Commit:

```bash
git add docs/compatibility.md scripts/capture_dwsim_reference.ps1 tests/golden
git commit -m "test: freeze initial DWSIM reference cases"
```

**Exit gate:** Golden cases are reproducible and include source revision, model settings, values, and units.

### Phase 1: Create the minimal Python package

**Files:**
- Create: `pyproject.toml`
- Create: `src/mesim/__init__.py`
- Create: `src/mesim/errors.py`
- Create: `scripts/validate.py`

**Steps:**

1. Add a Python 3.12 package with no runtime dependency.
2. Add domain errors: `ValidationError`, `MissingCompoundData`, `OutOfRangeError`, and `ConvergenceError`.
3. Make `scripts/validate.py --quiet` supported from the first commit.
4. Write a failing self-check proving the package imports and validation returns exit `0` with no enabled cases.
5. Implement only enough package metadata and validation plumbing to pass.
6. Run `python scripts/validate.py --quiet` and confirm exit `0`.
7. Commit `chore: create cross-platform calculation package`.

**Exit gate:** Imports and validation work on macOS arm64 and Linux in CI.

### Phase 2: Units of measure

**Files:**
- Create: `src/mesim/units.py`
- Create: `tests/test_units.py`

**Steps:**

1. Write failing checks for temperature, pressure, mass flow, molar flow, energy, power, enthalpy, density, viscosity, thermal conductivity, area, volume, length, velocity, and heat-transfer coefficient.
2. Add a fixed unit table and explicit `to_si(value, unit, dimension)` and `from_si(si_value, unit, dimension)` functions.
3. Reject unknown units, wrong dimensions, non-finite values, and temperatures below absolute zero.
4. Preserve the original value and unit in immutable `Quantity` records.
5. Run `python -m unittest tests.test_units -v`.
6. Commit `feat: add explicit unit conversion boundary`.

Phase 2B extends this boundary with the process dimensions required by the current internal-alpha roadmap. Its design and implementation steps are recorded in `docs/plans/2026-07-13-phase-2b-units-design.md` and `docs/plans/2026-07-13-phase-2b-units.md`.

**Exit gate:** Every accepted kernel and API boundary dimension needed by the implemented phases is covered; equation and round-trip checks pass without hidden rounding. Later models add exact unit aliases only with failing calculation-boundary tests. Full DWSIM display-unit parity is not required.

### Phase 3: Compound records and versioned data

**Files:**
- Create: `src/mesim/compounds.py`
- Create: `data/compounds/v1.json`
- Create: `data/interactions/pr-v1.json`
- Create: `tests/test_compounds.py`

**Steps:**

1. Write failing checks for schema validation, CAS uniqueness, immutable IDs, exact units, provenance, and missing required PR properties.
2. Add only methane, ethane, propane, n-butane, and n-pentane.
3. Load records with the standard `json` module and immutable dataclasses.
4. Add symmetric PR binary interaction lookup with explicit missing-value behavior; default zero is allowed only when the dataset explicitly says so.
5. Compare the five records against captured DWSIM reference JSON and the cited source.
6. Commit `feat: add versioned light-hydrocarbon data`.

**Exit gate:** Compound and interaction records are deterministic, source-backed, immutable, and sufficient for PR.

### Phase 4: Ideal properties and correlations

**Files:**
- Create: `src/mesim/thermo/ideal.py`
- Create: `tests/test_thermo.py`

**Steps:**

1. Write failing checks for ideal-gas enthalpy, entropy, heat capacity, density, and selected vapor-pressure correlations.
2. Implement correlations only for the five initial compounds and only within documented temperature ranges.
3. Integrate heat capacity analytically where the correlation permits; use SciPy quadrature only if a required correlation cannot be integrated safely.
4. Reject extrapolation unless the caller explicitly opts in and the result carries a warning.
5. Compare equation cases and DWSIM results.
6. Commit `feat: add ideal compound properties`.

**Exit gate:** Properties pass independent equation cases and range enforcement.

### Phase 5: Peng-Robinson EOS

**Files:**
- Create: `src/mesim/thermo/peng_robinson.py`
- Modify: `tests/test_thermo.py`

**Steps:**

1. Write failing pure-component checks for `a`, `b`, alpha, cubic roots, fugacity coefficient, density, and departure properties.
2. Implement pure-component PR.
3. Write failing mixture checks for mixing rules and binary interactions.
4. Implement mixture PR and explicit vapor/liquid root selection.
5. Add stability and invalid-state checks; never choose a root only by array position without phase criteria.
6. Add NumPy only if the standard-library cubic implementation fails the equation test matrix; document the reason in the commit.
7. Compare against DWSIM PR cases across vapor, liquid, near-critical, and two-root regions.
8. Commit `feat: add verified Peng-Robinson properties`.

**Exit gate:** Fugacity, compressibility, density, enthalpy and entropy meet T1 tolerances for the five-compound domain.

### Phase 6: Flash calculations

Detailed execution plan: [2026-07-13-phase-6-flash-calculations.md](2026-07-13-phase-6-flash-calculations.md).

**Files:**
- Create: `src/mesim/thermo/flash.py`
- Create: `tests/test_flash.py`

**Steps:**

1. Write failing Rachford-Rice cases for all-liquid, all-vapor and two-phase states.
2. Implement bounded solution with explicit residual and iteration limit.
3. Write failing TP flash cases using PR fugacity iteration and phase stability checks.
4. Implement TP flash with a documented fallback path.
5. Add isothermal bubble-pressure and dew-pressure calculations with caller-supplied bounds.
6. Add PH using a bracketed temperature solve and an explicit enthalpy reference; defer PS until the expander needs it and defer TV/PV until a UO needs them.
7. Verify mass balance, fugacity equality, energy consistency, determinism, and DWSIM parity.
8. Commit `feat: add PR vapor-liquid flash calculations`.

**Exit gate:** All initial flash modes converge or return structured failure; no unconverged state is accepted.

### Phase 7: Streams and the first unit-operation slice

**Files:**
- Create: `src/mesim/streams.py`
- Create: `src/mesim/unitops/basic.py`
- Create: `tests/test_unitops.py`

**Steps:**

1. Add validated immutable `StreamState` inputs and calculated `PhaseState` outputs. Preserve ordered compound IDs; never align compositions by array position without checking IDs.
2. Implement explicit `flash_stream()`; constructing or copying a stream performs no hidden flash.
3. Implement mixer and splitter first because they require no new thermodynamic mode. Mixer requires an explicit outlet pressure not above the lowest inlet pressure, sums component flows and enthalpy flow, then performs PH flash. Splitter preserves state and composition exactly and changes only outlet flows from explicit fractions. Commit each separately after component and total balance checks.
4. Add an energy-stream scalar in SI power with one sign convention: positive duty enters the material stream.
5. Implement one explicit specification mode per thermal operation: outlet temperature for heater/cooler. Duty is calculated output; reject a negative heater duty or positive cooler duty. Additional modes wait for a real case.
6. Implement valve as fixed outlet pressure below inlet pressure with PH flash. Implement equilibrium separator only as routing of an already converged flashed state; outlet phase flow equals feed flow times phase fraction.
7. For each operation, test zero flow, invalid pressure, missing or conflicting specification, phase change, component balance, and energy closure.
8. Compare each operation to one deterministic DWSIM golden case. Commit one operation at a time.

**Exit gate:** U0 operations pass conservation checks and DWSIM parity without a flowsheet framework.

### Phase 8: Flowsheet graph and steady-state execution

**Files:**
- Create: `src/mesim/flowsheet.py`
- Create: `tests/test_flowsheets.py`

**Steps:**

1. Write failing checks for duplicate tags, missing ports, incompatible connections, multiple writers, cycles, and deterministic topological order. All cycles are rejected until the recycle phase exists.
2. Represent the flowsheet with plain dictionaries and dataclasses; do not introduce a graph dependency. Define fixed named ports per operation and validate connections before calculation.
3. Treat unit operations as pure calculations from immutable inlet states and specifications. Execute acyclic flowsheets in deterministic topological order with tag ordering as the tie-breaker.
4. Calculate into a new immutable result map. Publish nothing unless every operation converges; a failed solve retains the last complete authoritative state unchanged.
5. Add a four-operation golden flowsheet and compare every stream.
6. Commit `feat: solve acyclic steady-state flowsheets`.

**Exit gate:** A feed-heater-valve-separator flowsheet matches DWSIM and never partially mutates on failure.

### Phase 9: HTTP API, containers, and internal alpha

**Files:**
- Create: `src/mesim/api.py`
- Create: `tests/test_api.py`
- Create: `Dockerfile`

**Steps:**

1. Expose only the verified compound lookup, flash, unit-operation calculation, flowsheet validation, and solve functions.
2. Use Pydantic only at the HTTP boundary. Translate once into immutable kernel records and preserve submitted values and units beside SI values.
3. Reject unknown models, compounds, units, tags, duplicate IDs, and non-finite numbers before calculation.
4. Return structured solver failure records without stack traces or partial authoritative outputs.
5. Version the API schema independently of model and data versions.
6. Keep the service stateless. Set request-size, iteration, and execution-time limits at the API boundary.
7. Record the internal authentication and network-exposure decision before deployment; do not expose an unauthenticated calculation service publicly.
8. Build and test Linux arm64 and amd64 images. The production image contains no DWSIM files or Windows runtime.
9. Run the same normalized request vectors on macOS arm64 and both Linux architectures.

**Exit gate:** The first production-usable slice is callable from the internal web app and returns equivalent deterministic JSON on macOS arm64 and Linux amd64/arm64.

### Phase 10: Pressure-changing and vessel operations

**Files:**
- Create: `src/mesim/unitops/pressure.py`
- Create: `src/mesim/unitops/separation.py`
- Modify: `tests/test_unitops.py`

**Order:** pump, compressor, expander, component separator.

Implement one mode first: specified outlet pressure and constant efficiency. Pump requires liquid inlet; compressor requires vapor inlet; expander adds PS flash before implementation. Reject efficiency outside `(0, 1]`, invalid phase, and impossible pressure direction. For each operation, write one normal and one boundary case, verify component/energy balance and DWSIM parity, then commit. Add performance curves only after a captured case requires them.

Do not add pass-through steady-state tank or vessel classes; the existing equilibrium separator already covers that calculation. Add vessel holdup only with dynamics or an approved residence-time/sizing requirement.

**Exit gate:** U1 matrix is green for the supported PR domain.

### Phase 11: Heat exchangers

**Files:**
- Create: `src/mesim/unitops/heat.py`
- Modify: `tests/test_unitops.py`

**Order:** specified-duty exchanger, specified-UA exchanger, effectiveness/NTU, pinch checks, shell-and-tube thermal rating.

Keep thermal rating and mechanical design separate. Mechanical code sizing, filters, and solids separation are outside this phase. Report crossed temperatures, invalid LMTD, phase-change segments, and non-convergence explicitly. Require segment energy closure when phase change occurs.

**Exit gate:** U2 energy closure and temperature profiles pass golden cases.

### Phase 12: Hydraulics and relief calculations

**Files:**
- Create: `src/mesim/unitops/hydraulics.py`
- Create: `tests/test_hydraulics.py`

**Order:** Darcy-Weisbach single phase, fittings, elevation, pipe thermal coupling, orifice, Lockhart-Martinelli, Beggs-Brill, then relief valve. Petalas-Aziz is added only when a required case falls outside the accepted Beggs-Brill domain.

Every correlation must declare units, regime, valid range, primary source, and behavior outside range. Preserve roughness, diameter, elevation, and segment data exactly as submitted. Two-phase hydraulics starts only after TP/PH flashes and phase properties are verified. Relief calculations require an explicit standard and edition, scenario basis, relieving conditions, and independent benchmark; DWSIM agreement alone is insufficient for safety acceptance.

**Exit gate:** Pressure-drop profiles close and each enabled correlation passes published examples plus DWSIM comparison.

### Phase 13: Expand thermodynamic coverage by approved demand

**Files:**
- Extend: `src/mesim/thermo/`
- Extend: `data/compounds/` and `data/interactions/`
- Extend: `tests/test_thermo.py` and `tests/test_flash.py`

This is a capability backlog, not an automatic sequence. Approve one model because an accepted flowsheet requires it, then implement it on one branch. Likely dependency order:

1. PR78 only if the current five-compound domain needs its high-acentric-factor correction; SRK only for a required comparison.
2. IAPWS-IF97 before verified steam/water equipment cases.
3. Wilson or NRTL before the first nonideal liquid VLE case.
4. UNIQUAC or UNIFAC only when required interaction/group data is source-backed.
5. LLE/VLLE stability and flash before any decanter or liquid-liquid column.
6. Electrolytes, sour water, seawater, petroleum, black oil, and solids each require a separate approved domain plan and dataset.

Do not create a generic property-package registry until the second model exists. At that point use a dictionary from stable model ID to constructor, not a plugin framework.

**Exit gate:** Each model has a declared compound/data domain and its own passing compatibility slice before registration.

### Phase 14: Reactions and reactors

**Files:**
- Create: `src/mesim/unitops/reactors.py`
- Create: `data/reactions/`
- Create: `tests/test_reactors.py`

**Order:** explicit elemental-composition data, reaction schema and element balance, conversion reactor, equilibrium reactor, Gibbs reactor, CSTR, PFR.

Do not parse molecular-formula strings as authoritative element data. Add versioned elemental composition, formation enthalpy, formation Gibbs energy, and reaction provenance before consuming them. Reject unbalanced stoichiometry at load time. Store original rate-expression units. Equilibrium/Gibbs reactors require verified chemical potentials for every participating phase. CSTR/PFR require an approved kinetic-rate schema and SciPy integration; they are not extensions of the conversion reactor.

**Exit gate:** U4 and U5 pass independent reaction cases and DWSIM flowsheets.

### Phase 15: Recycles, specifications and flowsheet convergence

**Files:**
- Modify: `src/mesim/flowsheet.py`
- Create: `tests/test_recycles.py`

**Order:** tear-stream recycle, energy recycle, adjust, specification.

Require explicit tear variables, scaling, convergence tolerances, and initial guesses. Start with direct substitution and bounded damping. Add Wegstein or Newton only when a captured case proves direct substitution insufficient. Record every residual vector, norm, damping value, and chosen algorithm. Never publish a partially converged flowsheet.

**Exit gate:** Recycle cases converge deterministically from documented initial guesses or fail with complete residual history.

### Phase 16: Columns

**Files:**
- Create: `src/mesim/unitops/columns.py`
- Create: `tests/test_columns.py`

**Order:** equilibrium-stage data model and balances, shortcut column, absorber, reboiled absorber, tridiagonal/bubble-point solver, Newton solver, sum-rates/inside-out only when a benchmark requires them.

Each stage must close component material, phase equilibrium, summation, and energy equations. Do not begin rigorous columns until the required VLE/LLE model and Phase 15 recycle convergence are verified. Store stage profiles and residual history on failure.

**Exit gate:** U6/U7 benchmark columns converge from documented initial estimates and match stage profiles within tolerance.

### Phase 17: Dynamics and controls

**Files:**
- Create: `src/mesim/unitops/dynamics.py`
- Create: `tests/test_dynamics.py`

Define holdup states, algebraic constraints, initialization, and event handling before choosing an integrator. Add tanks/vessels, heat exchangers, PID control, then dynamic reactors. Use SciPy `solve_ivp` only for systems reduced explicitly to ODE form. If a required model remains a DAE, evaluate an IDA-capable dependency instead of writing a solver or pretending `solve_ivp` supports it. Require fixed-step reproducibility mode for regression cases and separate adaptive integration for production accuracy.

**Exit gate:** Mass and energy accumulation close over every time step and dynamic golden trajectories pass.

### Phase 18: Optional specialty energy operations

**Files:**
- Create only after approval: `src/mesim/unitops/specialty.py`
- Create only after approval: `tests/test_specialty.py`

This phase is not scheduled by default. Add hydroelectric, wind, solar, fuel-cell, or electrolyzer models only after an internal use case supplies requirements and source data. Treat manufacturer/model constants as versioned inputs, not hidden constants.

**Exit gate:** U10 cases pass published model examples and DWSIM comparisons.

### Phase 19: Performance and optional native acceleration

**Files:**
- Create only after profiling: `benchmarks/`
- Create only after approval: `native/`

1. Define measured latency and throughput targets from the deployed internal API, then benchmark the operations that miss them.
2. Optimize Python data movement and algorithms first.
3. Use NumPy/SciPy vectorization second.
4. Add C/C++ only if an approved benchmark remains below the target throughput.
5. Keep the Python implementation as the reference and run identical golden vectors against both implementations.
6. Fall back to Python when the native extension is unavailable.

**Exit gate:** Native code is optional, measured, cross-platform, and numerically equivalent.

### Phase 20: Compatibility and release process

**Files:**
- Modify: `docs/compatibility.md`
- Create: `docs/model-limitations.md`
- Create: `docs/data-provenance.md`
- Create: `CHANGELOG.md`

1. Mark every matrix item `unsupported`, `partial`, `verified`, or `verified-with-difference`.
2. Publish validity ranges, data sources, tolerances, unsupported DWSIM modes, and known numerical differences.
3. Run the Python calculation suite against committed golden references on macOS arm64 and Linux amd64/arm64; DWSIM capture itself remains Windows-only.
4. Version compound data, model implementations, API schema, and solver settings independently.
5. Invoke this release process first after the Phase 9 internal-alpha gate and tag `0.1.0`. Publish later verified capability slices as minor releases; do not wait for universal DWSIM parity.

**Exit gate:** Release scope is explicit, reproducible, auditable, and deployable without Windows.

## 7. First production-usable milestone

Stop after Phase 9 and ship an internal alpha when all of these are true:

- Five light hydrocarbons are source-backed and versioned.
- Units are explicit and round-trip safely.
- Ideal properties and Peng-Robinson pass equation and DWSIM checks.
- TP and PH flashes converge across the accepted domain.
- Mixer, splitter, heater, cooler, valve and separator conserve material and energy.
- One acyclic flowsheet solves deterministically on an M1 Mac and Linux amd64/arm64.
- The internal web app can call the versioned API from a Linux arm64/amd64 container.
- Unsupported models fail explicitly instead of returning approximate values.

This milestone proves the architecture. Do not begin columns, electrolytes, dynamics, or C++ before it passes.

## 8. Per-task execution loop

Use this loop for every equation, model, flash mode, and unit operation:

1. Add one failing equation or golden check.
2. Run the exact test and confirm the expected failure.
3. Implement the smallest valid calculation.
4. Run the targeted test.
5. Run `python scripts/validate.py --quiet`.
6. Update `docs/compatibility.md` only after both pass.
7. Commit one coherent model or operation.

Do not batch several unverified models into one commit.

## 9. Decisions intentionally deferred

- Database: JSON is sufficient until concurrent compound editing exists.
- Plugin system: unnecessary until external teams need third-party models.
- Distributed solving: unnecessary until one process misses measured demand.
- C/C++ kernel: unnecessary until profiling proves Python/NumPy/SciPy insufficient.
- Exact DWSIM file import: build after the Python object model stabilizes; initial golden cases use normalized JSON.
- Web flowsheet editor: separate product work after the API is stable.

## 10. Current implementation sequence

1. Finish Phase 5 DWSIM parity and merge the PR implementation candidate.
2. Execute the detailed Phase 6 plan through TP, bubble/dew pressure, PH, and DWSIM gates.
3. Implement Phase 7 streams and unit operations one commit at a time.
4. Implement Phase 8 acyclic atomic flowsheet execution.
5. Implement Phase 9 API and multi-architecture containers, then release the internal alpha.
6. Stop for measured user feedback before approving any Phase 10+ capability.

Do not implement columns, reactions, specialty thermodynamics, dynamics, or native code merely to advance the matrix.
