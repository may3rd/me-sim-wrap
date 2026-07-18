# Compatibility and reference-capture record

## Scope

The rewrite uses the vendored DWSIM source under `dwsim-windows/` as a reference oracle only. The new calculation kernel must not import or load DWSIM assemblies.

DWSIM is licensed under GPLv3. This workspace is currently for internal use. Any future distribution of code derived from DWSIM must be reviewed for GPLv3 compliance before release.

## Release support matrix

Status terms are normative: `verified` means equation and executable golden gates pass in the stated domain; `verified-with-difference` means the documented numerical or model distinction is intentional; `partial` means only the named slice is supported; and `unsupported` means callers must not expect a result.

| Capability | Status | Release boundary |
|---|---|---|
| T0 ideal-gas and pure-component correlations | verified | Versioned catalog compounds within each saved correlation range |
| T1 Peng-Robinson pure and mixture properties | verified-with-difference | Classic PR; unmodified and explicit Peneloux density paths remain distinct |
| T2 TP, PH, bubble, and dew flash | verified-with-difference | Executable DWSIM parity gates pass; phase labeling, reference convergence, and caloric-model differences are documented |
| T3 caloric and transport extensions | partial | Named vapor/liquid correlations and captured mixtures only |
| T4 activity-coefficient liquid VLE | partial | Saved-source acetone/methanol NRTL activities and modified-Raoult bubble/dew pressures only |
| U0 streams, basic operations, and acyclic flowsheets | verified-with-difference | Phase-split flow uses the documented reference-flash tolerance |
| U1 pressure-changing and component-separation operations | verified | Captured PR pump, compressor, expander, and component-separator modes |
| U2 heat exchangers | partial | Duty, UA, efficiency, terminal pinch, and one vapor shell-and-tube rating slice |
| U3 hydraulics and relief utilities | partial | Captured pipe/orifice/correlation equations; not a relief-system design tool |
| U4 equilibrium and conversion reactors | partial | Named single/multiple reaction and vapor Gibbs slices |
| U5 kinetic reactors | partial | One ideally mixed CSTR and one supplied-profile PFR reaction family |
| U6/U7 columns | partial | Shortcut and fixed-thermodynamic profile solvers; no live fully energy-coupled rigorous solve |
| U8 recycle and logical blocks | verified | Direct-substitution material/energy recycle, scalar adjust, and bounded expression specification |
| U9 dynamics and controls | partial | Explicit holdup/HX/CSTR primitives and one DWSIM tank/PID trajectory; no general DAE solver |
| U10 specialty energy | partial | Source-equation solar, wind, and hydro gates only |
| Electrolytes, solids, petroleum assays, CAPE-OPEN, fuel cells, and electrolyzers | unsupported | No Python calculation path |

The detailed evidence and tolerances below define each row. [Model limitations](model-limitations.md) are part of this matrix, not optional guidance.

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

`tests/golden/u2-heat-duty-pr-eos.json` is the fixed-duty heat-exchanger reference: methane/ethane hot and cold feeds, each 2 mol/s at 500 kPa, with 400 K and 300 K inlet temperatures and a 2 kW heat transfer. Positive duty transfers energy from the hot inlet to the cold inlet at unchanged pressure. It has no external energy stream. Both outlet temperatures and enthalpies match within `1e-4` relative. The earlier `u1-heat-duty-pr-eos` pair remains an immutable audit artifact superseded by this explicitly configured U2 reference.

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

`tests/golden/u3-pipe-thermal-estimated-htc-global-irradiated-liquid-pr-eos.json` and its normalized repeat exercise DWSIM's `UseGlobalSolarRadiation=true` path with the saved flowsheet weather value of 1 kWh/m2 and absorption efficiency 0.001. `pipe_solar_irradiation_source` selects that global value instead of the saved 0.01 kWh/m2 local field, producing the same independently checkable 2.789834 kW absorbed rate per active segment. Both supplied and automatically recalculated liquid-property profiles pass the captured `3e-3` segment-duty, `1e-4` outlet-temperature, and `2e-3` or 35 W aggregate-energy gates. DWSIM 9.0.5's legacy pipe editor does not write the global checkbox to the profile; the reference sets the saved flag explicitly, reopens the case, verifies the checked state, and solves in DWSIM. Automatic internal liquid-property updates now cover every captured estimated-air, external-water, and soil branch. External-water properties, burial inputs, and irradiated inlet volumetric flow remain explicit; phase change in estimated-HTC profiles remains unsupported.

`tests/golden/u3-pipe-thermal-phase-change-pr-eos.json` and its normalized repeat cover a defined-heat liquid-to-two-phase PR pipe. The 2 MW saved specification advances four active 400 kW rows; DWSIM's stream energy is 1.5995935 MW after pressure coupling. The reference remains liquid for the first two active outlets and becomes two-phase for the final two, ending at 339.449473 K, 391489.647 Pa, and 0.0712191 molar vapor fraction. `pipe_defined_heat_pr_profile` now preserves each PH flash phase, vapor fraction, and molar enthalpy instead of rejecting the first non-liquid row. Segment temperatures match within `3e-6` relative, the outlet phase split uses a documented `5e-4` reference-flash tolerance, total molar flow remains exact, and the independent PR stream-energy closure matches within `5e-5` relative. Pressure endpoints are still explicit inputs, so fully coupled thermal/two-phase hydraulic iteration remains outside this gate.

Windows setup and capture steps are in [phase-5-dwsim-parity.md](phase-5-dwsim-parity.md).

## Phase 13 T4 activity-model status

`tests/golden/t4-nrtl-acetone-methanol-vle.json` and its normalized repeat isolate the `HP Azeotrope` stream from DWSIM 10's accepted pressure-swing acetone-column case. The capture enables bubble/dew calculation in memory, leaves the saved column unchanged, and returns without solve, object, or property-read errors. A first candidate using `HP Feed` was rejected because DWSIM's near-azeotropic result placed bubble pressure 17 Pa below dew pressure. The accepted stream is 69.518685 mol% acetone and 30.481315 mol% methanol at 388.288289 K; its 607.921595 kPa bubble pressure remains above its 596.554401 kPa dew pressure.

`data/correlations/nrtl-acetone-methanol-v1.json` preserves the official sample's embedded COCO vapor-pressure records and both directed ChemSep NRTL parameter records in their saved cal/mol basis. The loader converts those parameters explicitly to J/mol and rejects missing interactions, invalid units, duplicate directed pairs, invalid correlation ranges, and non-UTC provenance. `nrtl_activity_coefficients` follows DWSIM's saved-source NRTL equation, while `nrtl_bubble_pressure` and `nrtl_dew_pressure` close the scoped modified-Raoult equations with normalized incipient-phase compositions. Equation vectors match at floating-point precision; captured activities and both DWSIM envelope pressures match within `2.5e-3` and `2e-3` relative respectively. The small reference difference is retained because the saved stream activities do not exactly equal a fresh evaluation of the saved interaction records. This is a bounded binary liquid-VLE dependency for future column work, not a general NRTL TP/PH flash, excess-enthalpy model, LLE/VLLE stability solver, or live rigorous-column thermodynamic loop.

## Phase 14 U4 reaction status

`tests/golden/u4-conversion-reactor-isomerization.json` and its normalized repeat freeze DWSIM 9.0.5's official Chao-Seader isomerization test after a clean executable solve. The active conversion reaction is N-butane to Isobutane at 33% base-reactant conversion with a saved -9.2 MJ/kmol reaction heat. `data/reactions/v1.json` stores explicit C4H10 element counts plus ideal-gas formation enthalpy, Gibbs energy, and entropy on one J/kmol basis; the loader rejects unbalanced stoichiometry and reaction heats inconsistent with those formation records. `conversion_reactor` reproduces the component material balance before product flashing. DWSIM writes separate vapor and liquid reactor outlets, so their aggregated trace-component flows use a documented `2e-7` reference-flash tolerance while aggregate total molar flow retains a `1e-9` gate. Reactor phase splitting, Chao-Seader thermodynamics, and adiabatic outlet-temperature parity remain separate unsupported calculations.

`tests/golden/u4-equilibrium-reactor-steam-reforming-pr-eos.json` and its normalized repeat isolate the equilibrium-reactor branch of DWSIM 9.0.5's official Gibbs/equilibrium sample at 1000 K and 101325 Pa. The feed is 2 mol/s methane plus 3 mol/s water; the two vapor reactions form CO2/H2 and CO/H2 on a PR fugacity basis. The kernel evaluates DWSIM's Gibbs equilibrium constants from explicit formation H/G data and its midpoint ideal-Cp integrals, solves both extents simultaneously with a damped Newton method, and closes its log chemical-potential residual below `1e-9`. All outlet component flows and reactant conversions match within `2e-5` relative, both extents within `2e-6`, and the 400.515144 kW isothermal duty within `2e-6`. This gate covers isothermal vapor reactions with Gibbs-derived equilibrium constants; liquid reactions, reaction approaches, and adiabatic temperature iteration remain unsupported in `equilibrium_reactor`.

`tests/golden/u4-gibbs-reactor-steam-reforming-pr-eos.json` and its normalized repeat isolate the Gibbs-reactor branch of that sample. `gibbs_reactor` uses log component flows and elemental Lagrange multipliers to minimize the same PR vapor chemical potentials without requiring a declared reaction set. The kernel closes C/H/O balances below `1e-11 kmol/s` and the dimensionless stationarity residual below `1e-9`; the minimized Gibbs objective decreases from inlet to outlet. DWSIM's saved IPOPT penalty solution leaves an oxygen-atom residual of about `8.4e-7 kmol/s`, so the physically constrained kernel intentionally does not reproduce that imbalance: component flows and conversions use a documented `1.2e-3` reference tolerance, and the 400.477870 kW DWSIM duty uses `2e-4`. Multiphase Gibbs minimization, solids, adiabatic temperature iteration, and alternate optimizer parity remain unsupported.

## Phase 14 U5 kinetic-reactor status

`tests/golden/u5-cstr-ethylene-glycol-raoult.json` and its normalized repeat isolate a 1 m3 isothermal CSTR built from DWSIM 9.0.5's ethylene-glycol kinetic sample. The liquid feed is 20 mol% ethylene oxide and 80 mol% water at 328.15 K and 500 kPa. `data/reactions/v1.json` preserves the saved molar-concentration basis, `kmol/m3` concentration unit, `kmol/[m3.h]` rate unit, forward Arrhenius factor of 0.005, zero activation energy in J/mol, and explicit component orders. `continuous_stirred_tank_reactor` uses SciPy bounded least squares to close the ideally mixed material-rate equation with the supplied outlet liquid volumetric flow; its rate, extent, and reference reaction heat match within `3e-5` relative. DWSIM's saved outlet retains about a `1e-3` relative material residual on the very small ethylene-glycol product flow, so stream flow and conversion parity use that documented tolerance while the kernel closes its own material-rate residual below `1e-15 kmol/s`.

`tests/golden/u5-cstr-methanol-carbonylation-uniquac.json` and its normalized repeat extend that gate to the saved 100 m3 adiabatic methanol/carbon-monoxide CSTR. The reaction name was changed from `Main Reaction_kinetics` to `Main Reaction kinetics` before capture because DWSIM 9.0.5 incorrectly parses every CSTR property containing an underscore as an indexed built-in property; the original capture was rejected after its rate, heat, and extent reads raised index errors. The repaired captures solve without object or property-read errors. Reaction data preserve the original `mol/m3` concentration and `mol/[m3.s]` rate bases, 3.5e6 Arrhenius factor, 83.68 kJ/mol activation energy, and 1.0/0.5 methanol/carbon-monoxide orders. The kernel converts both amount prefixes explicitly to its kmol/m3/s basis and matches every nonzero outlet component plus methanol conversion within `2.1e-4` relative. DWSIM's reported 0.902165 mol/s reaction extent lags the approximately 0.903494 mol/s final stream change, so its internally related extent, rate, and heat fields use a separate `1.7e-3` reference tolerance while the Python material-rate residual remains below `1e-12 kmol/s`. UNIQUAC phase-property prediction, automatic adiabatic temperature coupling, and multiple simultaneous kinetic reactions remain outside this fixed-state gate.

`tests/golden/u5-pfr-ethylene-glycol-raoult.json` and its normalized repeat cover the official 1 m3 adiabatic PFR sample with the same reaction. `plug_flow_reactor` uses SciPy `solve_ivp` over reactor volume and linearly interpolates the supplied inlet/outlet temperature and volumetric-flow endpoints. Component flows, conversion, integrated extent, average rate, and the 82.956 W reference reaction heat match within `2e-4` relative while the kernel closes stoichiometric material balance below `1e-15 kmol/s`. DWSIM displays the saved PFR extent and rate properties with mol-based unit labels, but their numeric values close the captured stream and heat balances only on a kmol basis; the gate follows that coherent basis explicitly. Automatic flash-property profiles, full adiabatic energy coupling, pressure-drop integration, multiple simultaneous reactions, and alternate PFR discretizations remain unsupported.

## Phase 15 recycle and logical-block status

`tests/golden/u8-material-recycle-cavett-pr-eos.json` and its normalized repeat isolate `REC-000` from DWSIM 9.0.5's official Cavett sample. The saved block uses acceleration method `None`, a 100-iteration limit, and explicit 1 K, 0.1 Pa, and 0.01 kg/s convergence tolerances. After the executable solve, inlet and outlet tear streams agree exactly in temperature, pressure, mass flow, specific enthalpy, and all 15 component molar flows. `solve_recycle` implements the approved first algorithm: damped direct substitution over explicit tear variables with caller-supplied scales, physical tolerances, damping, maximum iterations, and initial guess. Every evaluation records the guess, calculated vector, residual, scaled norm, damping, and iteration number; exhaustion raises `RecycleConvergenceError` carrying the complete immutable history and returns no partial flowsheet result. A separate contraction-vector test exercises deterministic multi-iteration convergence because the saved DWSIM case begins from its already converged tear state.

`tests/golden/u8-energy-recycle-turboexpander.json` and its normalized repeat isolate `EREC-116` from DWSIM's official natural-gas turbo-expansion sample. The controlled saved copy uses acceleration method `None`, a 100-iteration limit, and a 0.1 kW power tolerance. Both executable captures solve without object or property-read errors; the expander-minus-compressor power difference equals the reported recycle residual and is below tolerance, while the inlet and outlet energy streams carry the same accepted tear value. `solve_energy_recycle` fixes the scalar API to a coherent watt basis and delegates to the same bounded, fully recorded direct-substitution kernel.

`tests/golden/u8-adjust-biodiesel-nrtl.json` and its normalized repeat isolate `ADJ-000` from DWSIM's official NRTL Biodiesel Production sample. The adjust manipulates the `EtOH` feed mass flow in kg/s to drive `PROP_MS_104/Ethanol_BD` on the `Etanol` product to 1.6666666667 mol/s with a 0.0001 mol/s tolerance. Both executable captures solve without object or property-read errors and close the target to about `1.8e-8 mol/s`. `solve_adjust` implements bounded scalar finite-difference Newton steps with explicit initial guess, physical bounds, controlled-variable scale, tolerance, probe step, damping, and iteration limit. Each primary iteration records the manipulated and controlled values, target residual, scaled norm, derivative, applied step, and damping; singular derivatives, pinned bounds, and iteration exhaustion raise `AdjustConvergenceError` with the complete immutable history and return no partial result.

`tests/golden/u8-specification-petroleum-pr-eos.json` and its normalized repeat extend DWSIM's official PR petroleum-distillation sample with isolated energy-stream specification vectors. A source of 1234.56789 kW drives an affine `1.5 * X + 10.0` target to 1861.851835 kW and a bounded `-X` target to its saved -1000 kW lower limit. Both executable captures solve with no object or property-read errors. `apply_specification` evaluates a bounded, case-insensitive subset of DWSIM's `X` source, `Y` current-target, arithmetic, and `System.Math` expression model without exposing Python evaluation. It requires both bounds or neither, rejects unsupported syntax and non-finite states as validation errors, and records both unconstrained and applied results. Source and target values must use one caller-selected unit basis, matching DWSIM's selected-unit-system behavior. Accelerated recycle algorithms remain unsupported.

## Phase 16 column status

`equilibrium_stage_residuals` provides the first rigorous-column primitive on one coherent kmol/s and J/kmol basis. It preserves independent component-material, phase-equilibrium, liquid/vapor summation, and energy residuals for arbitrary trial states, and validates closure against caller-selected physical tolerances.

`tests/golden/u6-shortcut-column-methanol-pr-eos.json` and its normalized repeat freeze `SC-078` from DWSIM 9.0.5's official PR methanol-synthesis sample. The saved total-condenser case uses Methanol as the light key, Water as the heavy key, 0.1 light-key mole fraction in bottoms, 0.001 heavy-key mole fraction in distillate, reflux ratio 1.5, and 101325 Pa at both ends. `shortcut_column` reproduces DWSIM's Fenske non-key distribution, its ordinary Underwood-root branch, Gilliland stage count, feed-stage estimate, product rates and compositions, and four internal molar flows. The captured source stream's eight trial component fractions sum to 0.9999404093; DWSIM passes that vector directly to the shortcut routine, so the parity API deliberately does not normalize it. DWSIM also accepts its first product-rate iteration at a `1e-4` relative criterion while retaining the preceding bottoms rate; the kernel preserves that executable behavior rather than silently repairing the reference balance. Minimum reflux, actual stages, and feed stage match within `2e-6`, `3e-6`, and `3e-6` absolute respectively; product and internal flow values match at floating-point precision. The captured condenser and reboiler duties are gated against their attached energy streams, but duty and endpoint-temperature calculation are not yet implemented. Distributed non-key Underwood roots and fully energy-coupled rigorous column solves remain unsupported.

`tests/golden/u6-absorber-simple-pr-eos.json` and its normalized repeat freeze DWSIM's official six-component PR absorber. The saved 12-stage `Absorber` uses Burningham-Otto Sum-Rates with 100 iterations and `1e-10` internal and external tolerances. Both captures solve without object or property-read errors, remove more than 99.9999999% of feed propane, and close every component through `column_balance_residuals` within `2e-12 kmol/s`. The old sample stores localized stream keys (`Metano`, `Etano`, `Propano`, `nOctano`, `nNonano`, `nDecano`) while current DWSIM catalog lookup uses the corresponding English records; the golden preserves both exactly.

`fixed_k_sum_rates_absorber` is the predictive material-loop gate for that reference. Given the saved stage equilibrium-ratio matrix and initial liquid/vapor profiles, it reproduces the final liquid and vapor flows within `5.2e-7 kmol/s` and `1.4e-7 kmol/s`, and the liquid and vapor mole fractions within `5e-7` and `5e-6` absolute. It uses a tridiagonal component-flow solve on one kmol/s basis, applies DWSIM's component floor after converting it from mol/s, and preserves immutable iteration history on convergence failure. Energy-coupled temperature and equilibrium-ratio iteration remain unsupported.

`tests/golden/u6-reboiled-absorber-acetone-nrtl.json` and its normalized repeat freeze the 20-stage high-pressure column from DWSIM's official NRTL extractive-distillation sample after converting it to a reboiled absorber. The controlled case retains its `0.5 mol/s` bottoms specification, Wang-Henke bubble-point solver, experimental internal estimates, 1000-iteration limit, and `0.001` internal/external tolerances. The saved mode is explicitly `ReboiledAbsorber=true`; the former liquid-distillate mapping is removed and `HP Azeotrope` is saved as `OverheadVapor`. This mapping is required for DWSIM's executable material-balance path and is gated from the saved XML rather than inferred from the UI checkbox. Both captures reload and solve with no warnings, solve errors, object errors, or property-read errors. Total molar flow closes exactly and the two component-flow residuals are below `1.5e-9 kmol/s`.

The optional `-CaptureColumnProfiles` reference-capture switch preserves DWSIM's full-precision public `Tf`, `Lf`, `Vf`, `xf`, `yf`, and `Kf` arrays plus liquid and vapor `Hlf`/`Hvf` stage enthalpies in kJ/kg; it is opt-in so existing golden schemas stay stable. `fixed_k_material_column` uses a frozen thermodynamic closure and independently solves every stage component balance plus active-vapor summation with positive log flows, reduced softmax liquid compositions, finite-difference Newton steps, and a residual-reducing line search. Explicit liquid and vapor product draws and fixed zero-vapor endpoints cover total condensers. The reboiled-absorber gate reaches a scaled residual below `1.4e-15`; liquid and vapor compositions match within `4e-14`, while internal flows match within `2.2e-7 kmol/s`, consistent with correction of DWSIM's saved `0.001` column tolerance. Singular Jacobians, failed line searches, and iteration exhaustion return no partial profile and preserve immutable accepted-iteration history.

`tests/golden/u7-distillation-acetone-nrtl.json` and its normalized repeat freeze the original 20-stage, 6-atm acetone column from DWSIM 10's official extractive-distillation sample. The saved case is a total condenser with reflux ratio 40, a 0.5 mol/s bottoms specification, one feed on stage 10, Wang-Henke bubble-point solution, experimental internal estimates, 1000 iterations, and `0.001` internal/external tolerances. Both captures solve without solve, object, or property-read errors. The Python mesh gate closes material, phase equilibrium, and both summations on every stage within that saved tolerance; its maximum scaled component residual is below `3.3e-4`. DWSIM's public `Kf` trails `yf/xf` by as much as `5.2e-4` near the reboiler, so the Newton acceptance gate uses the captured phase ratio and separately audits public-`Kf` equilibrium closure. Starting from the captured internal estimate, the total-condenser profile is already converged at the saved tolerance and retains an exact zero overhead-vapor flow. `column_profile_energy_residuals` converts each composition-weighted kJ/kg stage enthalpy to J/kmol on the same kmol/s flow basis; all 20 stage heat residuals are below 31 W, or `1.5e-5` of the approximately 2.07 MW duty scale. The independent external-stream column balance closes within `2e-5 W`. Temperature, pressure-profile, duty/specification, and live NRTL equilibrium-ratio iteration remain unsupported, so Phase 16 establishes executable fixed-thermodynamic-profile parity rather than claiming a fully energy-coupled predictive rigorous-column solver.

## Phase 17 dynamics and controls status

`HoldupState`, `HoldupRates`, and `advance_holdup` define component accumulation on one kmol basis and energy accumulation on a J/W basis, returning explicit per-step residuals. Algebraic constraints and consistent initialization are represented and validated separately. `fixed_step_explicit_euler` supplies the deterministic regression path with grid-aligned state events, while `adaptive_explicit_ode` uses SciPy `solve_ivp` only for models already reduced to explicit ODE form. It rejects an asserted DAE model with an IDA-capable-solver requirement instead of treating a general algebraic system as an ODE.

`tests/u9-water-tank-level-control.dwxmz` freezes DWSIM 10's official two-cubic-metre water-tank level-control sample. `tests/golden/u9-water-tank-level-control-dwsim.txt` and its byte-identical repeat contain all 121 records at five-second intervals from two clean dynamic runs. The predictive gate reproduces DWSIM's explicit tank inventory update, hydrostatic valve pressure, liquid Kv equation, reverse-acting PID, integral guard, and one-integration-step outlet-flow lag. Every recorded level matches within `1e-4 m`, every outlet-valve position within `0.013` percentage point, and the Python mass accumulation residual remains below `4e-13 kg` per step. The position tolerance includes DWSIM's results-table rounding to four decimals.

The same phase adds equal-and-opposite UA heat transfer for a two-holdup lumped heat exchanger and a dynamic CSTR step driven by explicit stoichiometry, reaction extents, and reaction enthalpies. Unit gates close their total energy and component balances at floating-point precision. Multiphase dynamic flashes, pressure-flow networks, equipment geometry beyond this lumped tank, adaptive event location, and a production DAE solver remain unsupported.

## Phase 18 specialty energy status

`tests/u10-solar-panel.dwxmz`, `tests/u10-wind-turbine.dwxmz`, and `tests/u10-hydroelectric-turbine.dwxmz` freeze DWSIM 10's official renewable-energy samples; the solar case was solved and saved again through the executable before the gate was added. `specialty.py` implements the corresponding vendored DWSIM source equations with all weather and thermodynamic properties passed explicitly. It contains no live weather lookup or hidden manufacturer curve.

The saved global solar irradiation of `0.810493563416576 kW/m2`, one square-metre panel, and 15% efficiency produce `0.121574034512486 kW`. The wind gate uses the saved `1.17454042760979 kg/m3` humid-air density and `5.6 km/h` wind speed with DWSIM's `8/27` coefficient, ten square metres of disk area, 50 turbines, and 80% downstream efficiency to reproduce `0.654969046052182 kW` theoretical and `0.523975236841745 kW` generated power. The hydro gate derives volumetric flow from the saved mass flow and liquid density, adds static and velocity head with DWSIM's `9.8 m/s2` gravity, and matches `7.60872279398433 kW`; its outlet enthalpy decrement closes shaft power at floating-point precision.

These gates are source-equation parity, not equipment design or energy-yield forecasts. Site resource distributions, wind power curves and cut-in/cut-out behavior, PV temperature/angle/degradation effects, hydraulic losses beyond the supplied efficiency, water electrolyzers, and PEM fuel cells remain unsupported. Manufacturer-specific constants must be versioned inputs before any of those models are added.

## Phase 19 performance status

`benchmarks/benchmark_release.py` measures representative TP, PH, dynamic, and HTTP flows with explicit p95 thresholds below half of the API's five-second calculation deadline. The initial run showed 2.0–2.3 second endpoint latency despite millisecond-scale TP and dynamics kernels because the API started a fresh Windows Python interpreter for every request. The API now reuses a two-worker spawned-process pool; a timeout terminates and discards the pool before returning HTTP 408, and the next calculation recreates it. This preserves the hard process boundary while removing per-request interpreter startup.

After that change, the Windows CPython 3.12.2 reference run measured p95 values of 14.82 ms for the TP endpoint, 154.66 ms for the PH-valve endpoint, and 332.93 ms for the complete U0 endpoint. All direct-kernel and endpoint gates pass, so Phase 19 adds no native extension. Deployed concurrent-load targets and macOS/Linux baselines still require measurements on their actual hosts; native acceleration remains conditional on one of those approved workloads missing its target after Python profiling.

## Phase 7 and 8 U0 status

The Python kernel implements immutable material and energy streams; mixer, splitter, heater, cooler, valve, and equilibrium separator calculations; and deterministic acyclic flowsheet execution. `tests/golden/u0-pr-c1-c5.json` is a repeatable DWSIM 9.0.4 PR flowsheet capture covering methane and n-pentane feeds, mixer, heater, valve, and separator. `tests/test_flowsheets.py` matches the mixer, heater, valve, and phase-product stream temperature, pressure, and total molar flow to `1e-5` relative. Phase-product flow uses `1e-4` relative because DWSIM's default flash convergence leaves a documented phase-split difference. The Python catalog key is `N-pentane`; it is mapped explicitly to DWSIM's captured `n-Pentane` record.

## Phase 6 flash parity status

`tests/golden/pr-flash.json` and `tests/golden/pr-flash-repeat.json` are repeatable DWSIM 9.0.5 captures of the methane/ethane PR flash domain. The capture has no solve, object, or property-read errors. Python and DWSIM agree on stable liquid and stable vapor states. DWSIM reports the one-root near-critical stream as vapor while Python reports `single`; Python does not manufacture a vapor split from one EOS root. The two-phase vapor fraction differs by `3.96e-5` relative: DWSIM's default `NestedLoops` solver exits when the vapor-fraction update is below `1e-6`, although the saved state has a fugacity residual of about `1.2e-4`; Python retains its `1e-8` fugacity-equilibrium requirement.

DWSIM's PR package is configured to use Lee-Kesler calorics, while this kernel uses ideal-gas heat-capacity correlations plus PR departure enthalpy. DWSIM-target PH temperatures differ by at most `2e-5` relative for the captured single-vapor and phase-crossing cases; Python PH energy closure remains `1e-6` relative to its supplied target.

The recapture enables DWSIM's `CalculateBubbleAndDewPoints` setting in memory without modifying the saved flowsheet. Both streams report the same nonzero envelope values at 180 K: approximately 2.162245 MPa bubble pressure and 0.262945 MPa dew pressure, with the required bubble-above-dew inequality. The Python PR pressure solvers match those properties within `2e-6` relative, while retaining their stricter internal fugacity residuals. This closes the T2 bubble/dew gate for the captured methane/ethane domain; it is not a general phase-envelope tracing capability.
