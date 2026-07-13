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

T1 equation checks cover vapor, liquid, single-root, three-root, and mixture states. Pure methane vapor and n-pentane liquid fugacity coefficients are checked against `tests/golden/u0-pr-c1-c5.json` from DWSIM 9.0.4. The T1 exit gate remains open until DWSIM references cover all five compounds, mixture fugacity, compressibility, density, departure enthalpy/entropy, near-critical states, and two-root states.

Windows setup and capture steps are in [phase-5-dwsim-parity.md](phase-5-dwsim-parity.md).
