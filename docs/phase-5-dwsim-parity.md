# Phase 5 DWSIM parity capture

Create one DWSIM 9.0.4 flowsheet and save it as `tests/pr-t1.dwxmz`.

## 1. Configure the flowsheet

Add the ChemSep records for these compounds and verify the constants before continuing:

| Compound | Tc (K) | Pc (Pa) | Acentric factor |
|---|---:|---:|---:|
| Methane | 190.56 | 4599000 | 0.011 |
| Ethane | 305.32 | 4872000 | 0.099 |
| Propane | 369.83 | 4248000 | 0.152 |
| N-butane | 425.12 | 3796000 | 0.199 |
| N-pentane | 469.7 | 3370000 | 0.251 |

Do not use the CoolProp n-butane or n-pentane records. Select `Peng-Robinson (PR)` as the property package and retain DWSIM's default flash algorithm.

## 2. Add reference streams

Add unconnected material streams with molar flow `1 kmol/h`. Use the tags, values, and SI units exactly as listed.

### Pure vapor

Each stream contains only the named compound.

| Tag | Compound | Temperature (K) | Pressure (Pa) |
|---|---|---:|---:|
| PR-V-METHANE | Methane | 228.672 | 459900 |
| PR-V-ETHANE | Ethane | 366.384 | 487200 |
| PR-V-PROPANE | Propane | 443.796 | 424800 |
| PR-V-NBUTANE | N-butane | 510.144 | 379600 |
| PR-V-NPENTANE | N-pentane | 563.64 | 337000 |

### Pure liquid

Each stream contains only the named compound.

| Tag | Compound | Temperature (K) | Pressure (Pa) |
|---|---|---:|---:|
| PR-L-METHANE | Methane | 133.392 | 2299500 |
| PR-L-ETHANE | Ethane | 213.724 | 2436000 |
| PR-L-PROPANE | Propane | 258.881 | 2124000 |
| PR-L-NBUTANE | N-butane | 297.584 | 1898000 |
| PR-L-NPENTANE | N-pentane | 328.79 | 1685000 |

### Special states

| Tag | Mole fractions | Temperature (K) | Pressure (Pa) |
|---|---|---:|---:|
| PR-3ROOT-METHANE | Methane 1.0 | 150 | 1000000 |
| PR-NC-METHANE | Methane 1.0 | 188.6544 | 4369050 |
| PR-MIX-ME-C2 | Methane 0.7; Ethane 0.3 | 250 | 5000000 |

Calculate the flowsheet. Stop if any stream reports a calculation error. Save the calculated file as `tests/pr-t1.dwxmz`.

## 3. Capture twice on Windows

```powershell
git fetch origin
git switch codex/phase-5-pr
git pull

$engine = 'C:\Users\26008353\Downloads\DWSIM_v904_win64_portable\Windows'

powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\scripts\capture_dwsim_reference.ps1 `
  -EngineBin $engine `
  -DwsimRevision '9.0.4' `
  -CasePath '.\tests\pr-t1.dwxmz' `
  -CaseId 'pr-t1' `
  -PropertyPackage 'Peng-Robinson (PR)' `
  -FlashAlgorithm 'DWSIM default' `
  -OutputPath '.\tests\golden\pr-t1.json'

powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\scripts\capture_dwsim_reference.ps1 `
  -EngineBin $engine `
  -DwsimRevision '9.0.4' `
  -CasePath '.\tests\pr-t1.dwxmz' `
  -CaseId 'pr-t1' `
  -PropertyPackage 'Peng-Robinson (PR)' `
  -FlashAlgorithm 'DWSIM default' `
  -OutputPath '.\tests\golden\pr-t1-repeat.json'
```

## 4. Verify and publish

```powershell
python .\scripts\validate.py --compare `
  .\tests\golden\pr-t1.json `
  .\tests\golden\pr-t1-repeat.json
```

Expected output:

```text
validation: normalized cases match
```

Do not commit if validation reports a mismatch or a property read error.

```powershell
git add tests/pr-t1.dwxmz tests/golden/pr-t1.json tests/golden/pr-t1-repeat.json
git commit -m "test: capture DWSIM Peng-Robinson references"
git push
```

