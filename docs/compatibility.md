# Compatibility and reference-capture record

## Scope

The rewrite uses the vendored DWSIM source under `dwsim-windows/` as a reference oracle only. The new calculation kernel must not import or load DWSIM assemblies.

DWSIM is licensed under GPLv3. This workspace is currently for internal use. Any future distribution of code derived from DWSIM must be reviewed for GPLv3 compliance before release.

## Source revision

The exact DWSIM source revision is not assumed from the directory name. Every capture must provide `-DwsimRevision` with the Git commit, release identifier, or an approved immutable snapshot identifier. The revision is stored in every golden case under `source.dwsim_revision`.

The current source areas used by Phase 0 are:

- `dwsim-windows/DWSIM.Automation/Automation.cs`
- `dwsim-windows/DWSIM.Interfaces/`
- `dwsim-windows/DWSIM.Thermodynamics/BaseClasses/ThermodynamicsBase.vb`
- `dwsim-windows/DWSIM.Thermodynamics/PropertyPackages/`
- `dwsim-windows/DWSIM.Thermodynamics/FlashAlgorithms/`
- `dwsim-windows/DWSIM.UnitOperations/UnitOperations/`
- `dwsim-windows/DWSIM.SharedClasses/UnitsOfMeasure/SystemsOfUnits.vb`

## Golden schema

Golden cases use `tests/golden/schema.json` with schema version `golden-case-1`. Each case records:

- immutable input values and units before calculation
- output values and units after calculation
- compound constants and provenance
- DWSIM revision and reported automation version
- property package and flash algorithm labels
- platform and architecture
- input flowsheet SHA-256 when a flowsheet file is used
- solver success, errors, object state, and property read failures

Numeric values include `value_text` in invariant culture so the source representation is retained alongside the JSON number.

## Unit conversion boundary

`src/mesim/units.py` uses an explicit allowlist. It does not parse arbitrary unit expressions and does not import DWSIM conversion code. Unknown symbols and physically incompatible dimensions are rejected.

The current canonical calculation bases are:

| Quantity family | Base unit |
|---|---|
| absolute temperature; temperature difference | K; K difference |
| pressure; pressure gradient | Pa; Pa/m |
| time; mass; amount | s; kg; kmol |
| mass, molar, and volumetric flow | kg/s; kmol/s; m3/s |
| energy; power; force | J; W; N |
| mass and molar enthalpy | J/kg; J/kmol |
| mass and molar heat capacity or entropy | J/kg/K; J/kmol/K |
| density; specific volume | kg/m3; m3/kg |
| dynamic viscosity; kinematic viscosity or diffusivity | Pa.s; m2/s |
| thermal conductivity; surface tension | W/m/K; N/m |
| length; area; volume; velocity; acceleration | m; m2; m3; m/s; m/s2 |
| mass flux; heat flux; heat-transfer coefficient | kg/m2/s; W/m2; W/m2/K |
| dimensionless values | 1 |

Gauge-pressure aliases are affine inputs referenced to exactly 101325 Pa and produce absolute Pa internally. Temperature differences never use absolute-temperature offsets. `Quantity` retains the submitted value and symbol separately from its derived calculation value.

Unit coverage grows only when an implemented model or approved API boundary needs it. Each addition requires an equation vector, round trip, and incompatible-dimension rejection. Matching every DWSIM display alias is not a compatibility target.

## Phase 0 capture command

Run this on x64 Windows with a built DWSIM output directory containing `DWSIM.Automation.dll` and its dependencies:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\capture_dwsim_reference.ps1 `
  -EngineBin 'C:\path\to\dwsim\bin\x64\Release' `
  -DwsimRevision 'REPLACE_WITH_COMMIT_OR_SNAPSHOT_ID' `
  -OutputPath '.\tests\golden\compound-catalog.json'
```

For a flowsheet case, add `-CasePath`, `-CaseId`, `-PropertyPackage`, and `-FlashAlgorithm`:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\capture_dwsim_reference.ps1 `
  -EngineBin 'C:\path\to\dwsim\bin\x64\Release' `
  -DwsimRevision 'REPLACE_WITH_COMMIT_OR_SNAPSHOT_ID' `
  -CasePath '.\PlatformFiles\Common\tests\basic\heating and cooling.dwxmz' `
  -CaseId 'basic-heating-cooling' `
  -PropertyPackage 'Peng-Robinson (PR)' `
  -FlashAlgorithm 'DWSIM default' `
  -OutputPath '.\tests\golden\basic-heating-cooling.json'
```

Capture the initial catalog plus these cases before Phase 1 starts:

1. single vapor
2. single liquid
3. two-phase flash
4. mixer
5. heater
6. valve
7. equilibrium separator

Use existing sample flowsheets when they isolate the operation. If a sample includes extra operations, record that in `source.notes` instead of silently treating it as an isolated test.

## Acceptance

Run each capture twice with the same binary and inputs. Compare the normalized JSON hashes. A mismatch blocks the phase until the nondeterminism is explained and either removed or recorded as an explicit limitation.

## Implemented thermodynamics

- T0: five-compound ideal-gas heat capacity, enthalpy, entropy, density, and pure vapor pressure.
- T1 implementation candidate: five-compound Peng-Robinson pure and mixture parameters, physical cubic roots, fugacity coefficients, density, departure enthalpy/entropy, and minimum-residual-Gibbs stable-root selection.

T1 equation checks cover vapor, liquid, single-root, three-root, near-critical, and mixture states. `tests/golden/pr-t1.json` and its normalized repeat capture DWSIM 9.0.4 references for all five compounds plus a methane/ethane mixture. Pure and mixture fugacity coefficients agree within `2e-6` relative. Reference-independent liquid-to-vapor enthalpy and entropy changes, including PR departure terms, agree within `2e-4` relative.

DWSIM material-stream density fields use `AUX_LIQDENS`/`AUX_VAPDENS`, while its PR compressibility path can apply configured Peneloux volume translation. These fields are recorded as `verified-with-difference`, not treated as equivalent to this implementation's unmodified PR cubic density and roots. Independent vapor and liquid equation vectors verify the Python PR compressibility and density calculations across all five compounds to `1e-8` relative. With that explicit model boundary, the classic-PR Phase 5 T1 gate is closed; PR78 remains demand-driven T3 scope.

T3A provides vapor dynamic viscosity and thermal conductivity for the five catalog compounds from ChemSep equation 102 curves. Conductivity is mole-averaged; viscosity is mole-averaged then corrected with DWSIM's Jossi-Stiel-Thodos dense-gas equation. `tests/golden/pr-t1.json` verifies those calculations against DWSIM PR values using the captured vapor density. The current unmodified PR density must not be substituted for a DWSIM volume-translated density in a parity claim; liquid transport and volume translation remain unsupported.

T3B adds an explicit Peneloux-translated vapor-density function using DWSIM's PR fallback coefficients for the same five compounds. It leaves classic PR states unchanged and matches the captured methane, ethane, and methane/ethane vapor densities within `1e-4` relative. Liquid transport remains unsupported.

T3C derives mixture Cp by a centered, fixed-pressure derivative of the existing converged TP-flash enthalpy. The methane/ethane PR reference at 250 K and 5 MPa agrees within `0.7%`; phase-change Cp and liquid transport remain unsupported.

`tests/golden/u1-pump-pr-eos.json` is the Phase 10 pump reference. Its DWSIM PR package sets liquid-density calculation to `EOS` and disables Peneloux volume translation; the default `Rackett_and_ExpData` density mode is not comparable because pump work depends on liquid molar volume. With the captured 75% efficiency, Python outlet temperature, enthalpy, and pump power match within `1e-4` relative. The earlier `u1-pump-pr*.json` captures remain immutable audit artifacts but are not pump-parity references.

`tests/golden/u1-compressor-pr-eos.json` is the PR compressor reference: methane/ethane vapor, 500 kPa to 1 MPa, DWSIM adiabatic path, and 75% adiabatic efficiency. The kernel uses a bounded PS flash for the isentropic outlet followed by PH flash for actual outlet enthalpy. It accepts converged `vapor` and one-root `single` inlets; a one-root state is not relabelled as vapor. Outlet temperature, enthalpy, and power match DWSIM within `1e-4` relative.

`tests/golden/u1-expander-pr-eos.json` is the matching adiabatic expander reference: methane/ethane, 1 MPa to 500 kPa, and 75% adiabatic efficiency. Temperature, enthalpy, and generated-power magnitude match within `1e-4` relative. `EnergyStream.duty_w` is negative because work leaves the material stream; DWSIM reports generated-power magnitude as positive.

`tests/golden/u1-component-separator-pr-eos.json` is the component-separator reference. Its saved specification is outlet 1, 10% n-pentane and 90% ethane of each component's inlet mass flow. Per-component mass and molar split fractions are numerically identical, so the kernel has one fraction-based mode. Both outlet flows and enthalpies match within `1e-4` relative. DWSIM reports `EnergyImb` as inlet minus outlets; `EnergyStream.duty_w` is its negation because positive duty enters material.

`tests/golden/u2-heat-duty-pr-eos.json` is the fixed-duty heat-exchanger reference: methane/ethane hot and cold feeds, each 2 mol/s at 500 kPa, with 400 K and 300 K inlet temperatures and a 2 kW heat transfer. Positive duty transfers energy from the hot inlet to the cold inlet at unchanged pressure. It has no external energy stream. Both outlet temperatures and enthalpies match within `1e-4` relative.

`tests/golden/u2-heat-ua-pr-eos.json` is the matching specified-UA reference: the same counter-current, no-loss streams with U = 25 W/m2/K and A = 1 m2. The kernel solves `Q = UA × LMTD` by bounded bisection and matches DWSIM's 1.955691 kW heat duty and both outlet temperatures within `1e-4` relative. Thermal-efficiency, pinch, phase-change, and general shell-and-tube modes remain unsupported.

`tests/golden/u2-heat-efficiency-pr-eos.json` covers DWSIM's `Specify Heat Transfer Efficiency` mode at 50%. The kernel calculates `Qmax` from both streams' inlet enthalpy changes at the opposite inlet temperature, uses the smaller value, then transfers the requested percentage. It matches DWSIM's 4.497993 kW duty and both 351.864 K outlet temperatures within `1e-4` relative. Pinch, phase-change, and general shell-and-tube modes remain unsupported.

`tests/golden/u2-heat-pinch-pr-eos.json` covers DWSIM's counter-current `Pinch Point` mode with MITA = 20 K and `PinchPointAtOutlets=false`. In the captured phase-free domain, the profile minimum is at a terminal; the kernel solves that minimum approach exactly. DWSIM's 25-segment, 0.01 kW profile search reports 7.086055 kW, while the kernel's exact result is accepted within 20 W; outlet temperatures match within `1e-4` relative. Phase-change profiles and the explicit outlet-pinch option remain unsupported.

`tests/golden/u2-shell-tube-rating-pr-eos.json` covers all-vapor, layout-0, counter-current shell-and-tube rating. It uses DWSIM's active steady-state baffle-derived shell-flow areas, Gnielinski tube coefficient, simplified Tinker shell correlation, external-area resistance stack, and corrected LMTD. Python matches the 6.461926 kW duty, both outlet temperatures, U, Reynolds numbers, and pressure drops within `1e-3` relative. Other tube layouts, liquid service, phase change, fouling, and pressure-coupled outlet flashes remain unsupported.

`tests/golden/u3-pipe-liquid-pr-eos.json` covers a five-increment, pure n-pentane PR pipe at 300 K. The kernel reproduces DWSIM's friction and static pressure-drop totals from each captured liquid segment state within `1e-5` relative, using its source friction-factor regimes: laminar below Re 2100, Churchill transition, and the explicit turbulent correlation above Re 4000. It reports velocity, Reynolds number, Darcy factor, friction drop, static drop, and total drop in SI. Fittings, thermal coupling, compressible flow, two-phase flow, stream-property calculation, relief sizing, and full pipe-operation parity remain unsupported.

Windows setup and capture steps are in [phase-5-dwsim-parity.md](phase-5-dwsim-parity.md).

## Phase 7 and 8 U0 status

The Python kernel implements immutable material and energy streams; mixer, splitter, heater, cooler, valve, and equilibrium separator calculations; and deterministic acyclic flowsheet execution. `tests/golden/u0-pr-c1-c5.json` is a repeatable DWSIM 9.0.4 PR flowsheet capture covering methane and n-pentane feeds, mixer, heater, valve, and separator. `tests/test_flowsheets.py` matches the mixer, heater, valve, and phase-product stream temperature, pressure, and total molar flow to `1e-5` relative. Phase-product flow uses `1e-4` relative because DWSIM's default flash convergence leaves a documented phase-split difference. The Python catalog key is `N-pentane`; it is mapped explicitly to DWSIM's captured `n-Pentane` record.

## Phase 6 flash parity status

`tests/golden/pr-flash.json` and `tests/golden/pr-flash-repeat.json` are repeatable DWSIM 9.0.4 captures of the methane/ethane PR flash domain. The capture has no property-read errors. Python and DWSIM agree on stable liquid and stable vapor states. DWSIM reports the one-root near-critical stream as vapor while Python reports `single`; Python does not manufacture a vapor split from one EOS root. The two-phase vapor fraction differs by `3.96e-5` relative: DWSIM's default `NestedLoops` solver exits when the vapor-fraction update is below `1e-6`, although the saved state has a fugacity residual of about `1.2e-4`; Python retains its `1e-8` fugacity-equilibrium requirement.

DWSIM's PR package is configured to use Lee-Kesler calorics, while this kernel uses ideal-gas heat-capacity correlations plus PR departure enthalpy. DWSIM-target PH temperatures differ by at most `2e-5` relative for the captured single-vapor and phase-crossing cases; Python PH energy closure remains `1e-6` relative to its supplied target.

Bubble and dew pressure remain **not DWSIM-parity verified**. The present `PR6-BUBBLE` and `PR6-DEW` streams merely store manually supplied pressure states and report `PROP_MS_126` and `PROP_MS_127` as zero because DWSIM's `CalculateBubbleAndDewPoints` setting is disabled. Do not mark T2 fully supported until a capture records nonzero DWSIM bubble/dew property values or a direct solver result.

Recapture on Windows with the capture switch below; it enables DWSIM's calculation setting in memory and does not modify the `.dwxmz` file:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\capture_dwsim_reference.ps1 `
  -EngineBin $engine `
  -DwsimRevision '9.0.4' `
  -CasePath '.\tests\pr-flash.dwxmz' `
  -CaseId 'pr-flash' `
  -PropertyPackage 'Peng-Robinson (PR)' `
  -FlashAlgorithm 'DWSIM default' `
  -CalculateBubbleAndDewPoints `
  -OutputPath '.\tests\golden\pr-flash.json'
```

Run it again with `pr-flash-repeat.json`, compare the two with `scripts/validate.py --compare`, and confirm nonzero `PROP_MS_126` and `PROP_MS_127` before committing the captures.
