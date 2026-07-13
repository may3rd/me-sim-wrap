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
| time; mass; amount | s; kg; mol |
| mass, molar, and volumetric flow | kg/s; mol/s; m3/s |
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
