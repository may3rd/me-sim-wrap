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

T3A provides vapor dynamic viscosity and thermal conductivity for the five catalog compounds from ChemSep equation 102 curves. Conductivity is mole-averaged; viscosity is mole-averaged then corrected with DWSIM's Jossi-Stiel-Thodos dense-gas equation. `tests/golden/pr-t1.json` verifies those calculations against DWSIM PR values using the captured vapor density. The current unmodified PR density must not be substituted for a DWSIM volume-translated density in a parity claim; volume translation remains an explicit model choice.

T3B adds an explicit Peneloux-translated vapor-density function using DWSIM's PR fallback coefficients for the same five compounds. It leaves classic PR states unchanged and matches the captured methane, ethane, and methane/ethane vapor densities within `1e-4` relative.

T3C derives mixture Cp by a centered, fixed-pressure derivative of the existing converged TP-flash enthalpy. The methane/ethane PR reference at 250 K and 5 MPa agrees within `0.7%`; phase-change Cp remains unsupported.

T3D extends the frozen five-compound transport catalog with ChemSep liquid-viscosity equation 101 and liquid-thermal-conductivity equation 16. `liquid_transport` follows the saved PR package configuration: mole-average pure viscosities and Li critical-volume-weighted conductivity mixing. The six liquid result rows in `tests/golden/u3-pipe-thermal-tabulated-htc-liquid-pr-eos.json` verify the 95.2381 mol% n-pentane / 4.7619 mol% ethane mixture from 300 to 301.501 K; both properties match DWSIM within `1e-12` relative. DWSIM extrapolates Ethane's conductivity curve above its saved 300 K maximum, while the Python API requires `allow_extrapolation=True` for those rows. The other three compounds have source-equation and catalog coverage but no separate executable liquid goldens.

`tests/golden/u1-pump-pr-eos.json` is the Phase 10 pump reference. Its DWSIM PR package sets liquid-density calculation to `EOS` and disables Peneloux volume translation; the default `Rackett_and_ExpData` density mode is not comparable because pump work depends on liquid molar volume. With the captured 75% efficiency, Python outlet temperature, enthalpy, and pump power match within `1e-4` relative. The earlier `u1-pump-pr*.json` captures remain immutable audit artifacts but are not pump-parity references.

`tests/golden/u1-compressor-pr-eos.json` is the PR compressor reference: methane/ethane vapor, 500 kPa to 1 MPa, DWSIM adiabatic path, and 75% adiabatic efficiency. The kernel uses a bounded PS flash for the isentropic outlet followed by PH flash for actual outlet enthalpy. It accepts converged `vapor` and one-root `single` inlets; a one-root state is not relabelled as vapor. Outlet temperature, enthalpy, and power match DWSIM within `1e-4` relative.

`tests/golden/u1-expander-pr-eos.json` is the matching adiabatic expander reference: methane/ethane, 1 MPa to 500 kPa, and 75% adiabatic efficiency. Temperature, enthalpy, and generated-power magnitude match within `1e-4` relative. `EnergyStream.duty_w` is negative because work leaves the material stream; DWSIM reports generated-power magnitude as positive.

`tests/golden/u1-component-separator-pr-eos.json` is the component-separator reference. Its saved specification is outlet 1, 10% n-pentane and 90% ethane of each component's inlet mass flow. Per-component mass and molar split fractions are numerically identical, so the kernel has one fraction-based mode. Both outlet flows and enthalpies match within `1e-4` relative. DWSIM reports `EnergyImb` as inlet minus outlets; `EnergyStream.duty_w` is its negation because positive duty enters material.

`tests/golden/u2-heat-duty-pr-eos.json` is the fixed-duty heat-exchanger reference: methane/ethane hot and cold feeds, each 2 mol/s at 500 kPa, with 400 K and 300 K inlet temperatures and a 2 kW heat transfer. Positive duty transfers energy from the hot inlet to the cold inlet at unchanged pressure. It has no external energy stream. Both outlet temperatures and enthalpies match within `1e-4` relative.

`tests/golden/u2-heat-ua-pr-eos.json` is the matching specified-UA reference: the same counter-current, no-loss streams with U = 25 W/m2/K and A = 1 m2. The kernel solves `Q = UA × LMTD` by bounded bisection and matches DWSIM's 1.955691 kW heat duty and both outlet temperatures within `1e-4` relative. Thermal-efficiency, pinch, phase-change, and general shell-and-tube modes remain unsupported.

`tests/golden/u2-heat-efficiency-pr-eos.json` covers DWSIM's `Specify Heat Transfer Efficiency` mode at 50%. The kernel calculates `Qmax` from both streams' inlet enthalpy changes at the opposite inlet temperature, uses the smaller value, then transfers the requested percentage. It matches DWSIM's 4.497993 kW duty and both 351.864 K outlet temperatures within `1e-4` relative. Pinch, phase-change, and general shell-and-tube modes remain unsupported.

`tests/golden/u2-heat-pinch-pr-eos.json` covers DWSIM's counter-current `Pinch Point` mode with MITA = 20 K and `PinchPointAtOutlets=false`. In the captured phase-free domain, the profile minimum is at a terminal; the kernel solves that minimum approach exactly. DWSIM's 25-segment, 0.01 kW profile search reports 7.086055 kW, while the kernel's exact result is accepted within 20 W; outlet temperatures match within `1e-4` relative. Phase-change profiles and the explicit outlet-pinch option remain unsupported.

`tests/golden/u2-shell-tube-rating-pr-eos.json` covers all-vapor, layout-0, counter-current shell-and-tube rating. It uses DWSIM's active steady-state baffle-derived shell-flow areas, Gnielinski tube coefficient, simplified Tinker shell correlation, external-area resistance stack, and corrected LMTD. Python matches the 6.461926 kW duty, both outlet temperatures, U, Reynolds numbers, and pressure drops within `1e-3` relative. Other tube layouts, liquid service, phase change, fouling, and pressure-coupled outlet flashes remain unsupported.

`tests/golden/u3-pipe-liquid-pr-eos.json` covers a five-increment, pure n-pentane PR pipe at 300 K. The kernel reproduces DWSIM's friction and static pressure drops from each captured liquid segment state within `1e-5` relative, using its source friction-factor regimes: laminar below Re 2100, Churchill transition, and the explicit turbulent correlation above Re 4000. `pipe_pressure_drop_profile` aggregates those supplied states and closes the captured friction total, static total, and outlet pressure to the same tolerance. `tests/golden/u3-pipe-two-phase-beggs-brill-pr-eos.json` and `tests/golden/u3-pipe-two-phase-lockhart-martinelli-pr-eos.json` verify all five saved segments for their selected DWSIM correlation: flow regime or holdup, friction, and hydrostatic drop. Their profile kernels reproduce every segment, the aggregate friction and static drops, and the product pressure within `1e-6` relative for Beggs-Brill and `1e-12` for Lockhart-Martinelli. It provides DWSIM's fixed-K fitting equation, its ISO-5167-style incompressible orifice calculation with corner, flange, and radius taps, its homogeneous Lockhart-Martinelli multiplier, its legacy-unit Beggs-Brill inclined-flow calculation when phase properties are supplied, and the DWSIM API RP 520 vapor-, liquid-, and two-phase-sizing utility equations. `tests/golden/u3-orifice-liquid-pr-eos.json` verifies the flange-tap n-pentane case: `2367.53195 Pa` plate drop and `1730.78032 Pa` recovered overall drop. `tests/golden/u3-psv-vapor-api520.json` records the methane PSV state and saved utility inputs; its confirmed DWSIM window result is `0.08013877 in²`, selecting `D / 0.11 in²`. `tests/golden/u3-psv-liquid-api520.json` records the liquid n-pentane state; the source equation returns `0.00480318 in²`, also selecting D. `tests/golden/u3-psv-two-phase-api520.json` records the 50/50 methane/n-pentane two-phase inlet and independent 90%-pressure state; the source equation returns `0.01905899 in²`, selecting D. The DWSIM utility does not expose a standard edition, so this is source-equation parity only, not relief-design or safety acceptance. Compressible flow, calculated phase properties, fitting catalogue data, relief discharge piping, and full pipe-operation parity remain unsupported.

`tests/golden/u3-pipe-thermal-liquid-pr-eos.json` and `tests/golden/u3-pipe-thermal-gradient-liquid-pr-eos.json`, each with a normalized repeat, are DWSIM 9.0.5 references for defined-HTC liquid-pipe thermal coupling. The captured feed is 95.2381 mol% n-pentane and 4.7619 mol% ethane at 300 K, 500 kPa, and 200 mol/s; the 100 m carbon-steel pipe has five saved increments, U = 25 W/m2/K, a 350 K base external temperature, and either zero or 0.1 K/m external-temperature gradient. `pipe_defined_htc_heat_transfer` uses the exact constant-property exponential form equivalent to the DWSIM LMTD equation when liquid heat capacity is supplied. `pipe_defined_htc_profile` advances that primitive across supplied segment lengths and heat capacities, while `pipe_defined_htc_gradient_profile` follows the vendored DWSIM source rule: evaluate ambient temperature from the segment-start cumulative distance and hold it constant through that increment. `liquid_pipe_supplied_state_profile` combines thermal and pressure profiles and enforces one mass-flow basis. The saved structures contain five hydraulic increments but advance the active thermal profile from result row 2 through row 5; the coupled gates preserve those distinct authoritative counts. In the gradient case, those starts are 20, 40, 60, and 80 m and the captured ambient temperatures are 352, 354, 356, and 358 K. Both cases match segment duties within `3e-3` relative, outlet pressure within `1e-5`, outlet temperature within `1e-4`, and energy within `2e-3`.

`tests/golden/u3-pipe-thermal-defined-heat-liquid-pr-eos.json` and its normalized repeat cover DWSIM's saved `Definir_Q` mode with a 10 kW specified load and five saved increments. The source assigns 2 kW to every result row. `pipe_defined_heat_pr_profile` applies that saved-increment divisor while advancing the active rows 2–5 with pressure-endpoint PH flashes, so its active four-row heat total is 8 kW and all outlet-row temperatures match within `2e-8` relative. The connected DWSIM energy stream is 7.593636 kW because its inlet-to-product balance also includes the omitted first hydraulic state; an independent PR inlet/product enthalpy difference matches that energy stream within `1e-5`. The active-row duty and full-stream energy are therefore gated separately rather than conflated.

`tests/golden/u3-pipe-thermal-tabulated-htc-liquid-pr-eos.json` and its normalized repeat cover defined-HTC mode with `UseUserDefinedU=true`. The saved table has distances 0, 50, and 100 m, external temperatures 350, 360, and 370 K, and U values 25, 30, and 35 W/m2/K. DWSIM linearly interpolates ambient temperature and applies a right-continuous step interpolation to U; at the active 20, 40, 60, and 80 m starts this produces 354, 358, 362, and 366 K with U = 25, 25, 30, and 30 W/m2/K. `pipe_tabulated_defined_htc_pr_profile` reproduces those rules and closes every increment with the captured pressure endpoint and one coherent kmol/s and J/kmol PR enthalpy balance. Segment duties match within `3e-3`, segment and product temperatures within `1e-4`, and aggregate energy within `2e-3` or 35 W absolute.

`tests/golden/u3-pipe-thermal-estimated-htc-liquid-pr-eos.json` and its normalized repeat extend that pipe to DWSIM's saved `Estimar_CGTC` mode at 350 K with 2 m/s external air. The authoritative saved switches include the internal Petukhov HTC, carbon-steel wall, and external Holman air HTC, with insulation and radiation disabled. `pipe_estimated_htc_air` reproduces the captured internal and external coefficients within `1e-12` relative, the wall coefficient within `5e-6`, and overall U within `1e-7`; `pipe_estimated_htc_air_profile` closes segment duties within `3e-3`, outlet temperature within `1e-4`, and energy within `2e-3`.

The same estimated-air reference now gates automatic liquid-property updates. `dwsim_pr_liquid_heat_capacity` implements the analytic `FluidProperties.PROPS.CpCvR` Peng-Robinson derivative used by DWSIM rather than approximating Cp from finite enthalpy differences. `pipe_liquid_pr_properties` combines that Cp with EOS density using DWSIM's source-level gas constant, the frozen experimental liquid-viscosity and Li conductivity correlations, and the saved pipe bore and molar-flow basis. Density, Cp, conductivity, viscosity, and velocity match all six captured rows within `1e-12` relative. `pipe_estimated_htc_air_pr_calculated_profile` recalculates those fields from each modeled segment inlet and retains the established `3e-3` segment-duty, `1e-4` outlet-temperature, and `2e-3` or 35 W aggregate-energy gates.

`tests/golden/u3-pipe-thermal-estimated-htc-gradient-liquid-pr-eos.json` and its normalized repeat add a 0.1 K/m ambient gradient to estimated-air-HTC mode. As in DWSIM's defined-HTC gradient, the active thermal rows start at 20, 40, 60, and 80 m and use 352, 354, 356, and 358 K external temperatures. `pipe_estimated_htc_air_pr_gradient_profile` applies those segment-start temperatures to both the Holman air properties and the pressure-endpoint PR enthalpy balance. `pipe_estimated_htc_air_pr_calculated_gradient_profile` additionally recalculates the five liquid correlation inputs from every modeled segment inlet while preserving the same start-distance rule. The liquid-domain bracket search stops at the first residual sign change rather than requiring the distant ambient endpoint to remain liquid. Internal and external HTC match within `1e-12`, wall HTC within `1e-5`, overall U within `1e-7`, segment duties within `3e-3`, outlet temperature within `1e-4`, and aggregate energy within `2e-3` or 35 W absolute.

`tests/golden/u3-pipe-thermal-estimated-htc-insulated-liquid-pr-eos.json` and its normalized repeat add the saved Fiberglass material record with 25 mm thickness and 0.035 W/m/K conductivity. Both supplied-property and automatically recalculated liquid-property profiles pass this resistance branch. Insulation and external-air coefficients match within `1e-12` relative, wall HTC within `5e-6`, overall U within `1e-7`, segment duties within `3e-3`, and outlet temperature within `1e-4`. At this low duty DWSIM's energy stream is 32 W below the sum of its own active result rows, so the aggregate energy gate uses a documented 35 W absolute tolerance. The deployed DWSIM 9.0.5 wall and insulation results use `k/(ln(Do/Di) Di)`-style resistances even though the vendored source snapshot contains an additional `/2`; executable golden behavior is authoritative for parity.

`tests/golden/u3-pipe-thermal-estimated-htc-water-liquid-pr-eos.json` and its normalized repeat switch the saved external medium to water at 310 K. The saved case is authoritative at 2 m/s external velocity, with insulation and radiation disabled. DWSIM 9.0.5 evaluates its IAPWS water properties at 310 K and 101325 Pa as 993.38984745012783 kg/m3 density, 0.00069354225499842137 Pa.s viscosity, 4178.7265325108107 J/kg/K heat capacity, and 0.62609581196362363 W/m/K conductivity; `pipe_estimated_htc_water` takes those external properties explicitly and reproduces the 5582.511740412966 W/m2/K external contribution within `1e-12`. Internal HTC also matches within `1e-12`, wall HTC within `1e-5`, and overall U within `5e-6`. `pipe_estimated_htc_water_pr_profile` closes supplied liquid properties at the captured pressure endpoints, while `pipe_estimated_htc_water_pr_calculated_profile` recalculates the internal liquid properties from each modeled segment inlet. Both use one coherent kmol/s and J/kmol PR enthalpy balance; segment duties match within `4e-3`, outlet temperature within `1e-4`, and aggregate energy within `2e-3`. Native IAPWS property calculation is outside this phase, so callers must still supply the external-water properties.

`tests/golden/u3-pipe-thermal-estimated-htc-dry-soil-liquid-pr-eos.json` and `tests/golden/u3-pipe-thermal-estimated-htc-moist-soil-liquid-pr-eos.json`, each with a normalized repeat, cover DWSIM's saved soil media at 310 K and 1 m initial burial depth. `dwsim_terrain_thermal_conductivity` freezes the vendored catalog mapping: gravel 1.1, stones 1.95, dry soil 0.5, and moist soil 2.2 W/m/K. `pipe_estimated_htc_soil` reproduces DWSIM's buried-cylinder resistance and the captured dry/moist external contributions of 2.4614250971428842 and 10.830270427428692 W/m2/K within `1e-12`; internal and wall coefficients retain their established `1e-12` and `5e-6` gates, and overall U matches within `1e-7`. `pipe_estimated_htc_soil_pr_profile` uses supplied liquid properties, while `pipe_estimated_htc_soil_pr_calculated_profile` recalculates them at every modeled segment inlet. Both share the pressure-endpoint PR enthalpy balance and match segment duties within `3e-3`, outlet temperature within `1e-4`, and aggregate energy within `2e-3` or 35 W absolute. Burial depth and conductivity remain explicit calculation inputs; gravel and stones have source-equation mapping coverage but no separate executable goldens.

`tests/golden/u3-pipe-thermal-estimated-htc-insulated-irradiated-liquid-pr-eos.json` and its normalized repeat enable local irradiation at 0.01 kWh/m2 with absorption efficiency 0.1. `pipe_absorbed_solar_radiation` follows DWSIM's steady source equation, including its use of the bare outer diameter for both residence volume and exposed area, yielding 2.789834 kW absorbed per active segment. `pipe_irradiated_heat_transfer` solves the supplied-property LMTD and radiation balance. The PR calculated-property driver adds the same absorbed rate inside each pressure-endpoint enthalpy residual and exposes both per-segment and aggregate radiation results. Both paths retain `3e-3` segment-duty, `1e-4` outlet-temperature, and 35 W aggregate-energy gates; inlet volumetric flow remains an explicit input because DWSIM's radiation source uses the inlet-stream residence time.

`tests/golden/u3-pipe-thermal-estimated-htc-global-irradiated-liquid-pr-eos.json` and its normalized repeat exercise DWSIM's `UseGlobalSolarRadiation=true` path with the saved flowsheet weather value of 1 kWh/m2 and absorption efficiency 0.001. `pipe_solar_irradiation_source` selects that global value instead of the saved 0.01 kWh/m2 local field, producing the same independently checkable 2.789834 kW absorbed rate per active segment. Both supplied and automatically recalculated liquid-property profiles pass the captured `3e-3` segment-duty, `1e-4` outlet-temperature, and `2e-3` or 35 W aggregate-energy gates. DWSIM 9.0.5's legacy pipe editor does not write the global checkbox to the profile; the reference sets the saved flag explicitly, reopens the case, verifies the checked state, and solves in DWSIM. Automatic internal liquid-property updates now cover every captured estimated-air, external-water, and soil branch. External-water properties, burial inputs, and irradiated inlet volumetric flow remain explicit, and phase change remains unsupported.

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
