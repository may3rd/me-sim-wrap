# Data provenance

All runtime property and reaction data is committed, versioned, and loaded without a DWSIM runtime dependency. Source timestamps are audit metadata; deterministic comparisons exclude capture timestamps.

The wheel installs these records under `share/mesim/data`; a source checkout uses the repository `data/` directory. Release smoke tests must import the installed wheel outside the checkout and confirm the same versioned records are found.

| Artifact | Schema | Primary source | Frozen revision | Verification |
|---|---|---|---|---|
| `data/compounds/v1.json` | `compound-data-1` | Vendored ChemSep `chemsep1.xml` | DWSIM 9.0.4 | 408 case-distinct fully supported records; loader schema, uniqueness, physical-field, and DWSIM captured-value tests |
| `data/interactions/pr-v1.json` | `pr-interactions-1` | `dwsim-windows/DWSIM.Thermodynamics/Assets/pr_ip.dat` | DWSIM 9.0.5.0 | 181 first-source pairs plus 16 accepted explicit-zero closures; symmetry, compound-key, and PR mixture gates |
| `data/interactions/srk-v1.json` | `srk-interactions-1` | `dwsim-windows/DWSIM.Thermodynamics/Assets/srk_ip.dat` | DWSIM 9.0.5.0 | 52 first-source pairs over the supported ChemSep domain; strict missing-pair and model-identity gates |
| `data/interactions/wilson-v1.json` | `dwsim-wilson-data-1` | Installed `DWSIM.Thermodynamics.wilson_bips.csv` resource plus the saved acetone/methanol Wilson case | DWSIM 9.0.5.0 | All 364 source pairs, embedded-resource hash, exact scoped 298.15 K molar volumes, and repeatable equation-level activity parity |
| `data/interactions/uniquac-v1.json` | `dwsim-uniquac-data-1` | Installed `DWSIM.Thermodynamics.uniquac.dat` resource plus the saved 1-propanol/water UNIQUAC case | DWSIM 9.0.5.0 | All 376 source pairs, embedded-resource hash, resolved binary `R`/`Q` basis, and repeatable equation-level activity parity |
| `data/interactions/unifac-v1.json` | `dwsim-unifac-data-1` | Installed original-UNIFAC group and interaction resources plus the saved 1-propanol/water case | DWSIM 9.0.5.0 | All 119 subgroups and 1,403 directed runtime interactions, both resource hashes, scoped surface-fraction vectors, and repeatable activity parity |
| `data/interactions/unifac-ll-v1.json` | `dwsim-unifac-data-1` | Installed original subgroup resource, UNIFAC-LL interaction resource, and saved 1-propanol/water case | DWSIM 9.0.5.0 | All 119 subgroups and 1,467 directed LLE interactions, both hashes, executable mixed-group normalization, and exact repeatable activity parity |
| `data/interactions/modfac-dortmund-v1.json` | `dwsim-modfac-data-1` | Installed Modified UNIFAC (Dortmund) group/interaction resources and saved 1-propanol/water case | DWSIM 9.0.5.0 | All 108 subgroups and 1,167 paired six-coefficient interactions, both hashes, scoped group vectors, and repeatable activity parity |
| `data/interactions/modfac-nist-v1.json` | `dwsim-modfac-data-1` | Installed Modified UNIFAC (NIST) group/interaction resources and saved 1-propanol/water case | DWSIM 9.0.5.0 | All 201 subgroups and 1,969 directed three-coefficient interactions, both hashes, scoped group vectors, and exact repeatable activity parity |
| `data/interactions/chao-seader-v1.json` | `dwsim-semi-empirical-data-1` | Installed Chao-Seader package over a saved methane/N-pentane case plus vendored package/model source | DWSIM 9.0.5.0 | Exact scoped pure constants, runtime and equation-source hashes, repeatable phase-fugacity vectors, and two-phase TP-flash parity |
| `data/interactions/grayson-streed-v1.json` | `dwsim-semi-empirical-data-1` | Installed Grayson-Streed package over a saved methane/N-pentane case plus vendored package/model source | DWSIM 9.0.5.0 | Exact scoped pure constants, separate equation-source hashes, repeatable phase-fugacity vectors, and two-phase TP-flash parity |
| `data/interactions/lkp-v1.json` | `dwsim-lkp-data-1` | Installed Lee-Kesler-Plöcker interaction resource, saved methane/N-pentane case, and vendored package/model source | DWSIM 9.0.5.0 | All 140 runtime interactions, both scoped critical-volume records, resource/runtime/source hashes, phase-fugacity parity, and two-phase TP-flash parity |
| `data/correlations/steam-tables-v1.json` | `dwsim-steam-tables-data-1` | DWSIM IAPWS-IF97 source and direct-class pure-Water golden case | DWSIM 9.0.5.0 | All ten region-4 coefficients, runtime/package/IAPWS hashes, exact executable fugacity parity, and pure TP phase classification |
| `data/thermo-packages/dwsim-v1.json` | `dwsim-property-package-catalog-1` | `dwsim-windows/DWSIM/Forms/MainWindow/FormMain.vb` (`AddPropPacks`) | DWSIM 9.0.5.0 | Exact 29-package built-in inventory; unique IDs/classes/names and runtime model-link gates |
| `data/correlations/ideal-v1.json` | `ideal-correlations-1` | Vendored ChemSep `chemsep1.xml` | DWSIM 9.0.4 | All 408 supported catalog compounds; range, midpoint, equation-vector, and representative DWSIM property gates |
| `data/correlations/transport-v1.json` | `transport-correlations-3` | Vendored ChemSep `chemsep1.xml` | DWSIM 9.0.4 | All 408 supported catalog compounds; positive midpoint checks and representative DWSIM phase transport parity |
| `data/correlations/saturated-liquid-v1.json` | `saturated-liquid-correlations-1` | Vendored ChemSep `chemsep1.xml` | DWSIM 9.0.4 | All 408 supported catalog compounds; positive midpoint checks and representative DWSIM density, heat-capacity, vaporization, and surface-tension parity |
| `data/correlations/chemsep-exclusions-v1.json` | `chemsep-exclusions-1` | Vendored ChemSep `chemsep1.xml` | DWSIM 9.0.4 | Audits all 23 excluded source records with exact missing fields and unsupported equation IDs |
| `data/correlations/nrtl-acetone-methanol-v1.json` | `nrtl-vle-data-2` | Saved COCO compound records, phase-caloric correlations, and ChemSep NRTL parameters inside the accepted pressure-swing column | DWSIM 10.0.9671.22371 case SHA-256 recorded in the file | Unit/schema rejection, equation vectors, deterministic NRTL envelope goldens, and all-stage caloric parity |
| `data/reactions/v1.json` | `reaction-data-1` | DWSIM isomerization, Gibbs/equilibrium, ethylene-glycol, and methanol-carbonylation reactor samples | DWSIM 9.0.5.0 | Element balance, formation-property consistency, original kinetic units, and reactor goldens |

## Golden references

- JSON references use `tests/golden/schema.json` (`golden-case-1`) and record DWSIM revision, executable version, platform, units, property-read errors, and input-flow SHA-256 where applicable.
- Captures span DWSIM 9.0.4, 9.0.5, and the committed DWSIM 10 sample revision identified inside each case. New captures must never replace an older oracle silently.
- `scripts/capture_dwsim_reference.ps1` is the Windows automation boundary. Invoke it with `-ExecutionPolicy Bypass`, an explicit engine directory, an explicit revision, and explicit case metadata.
- `tests/golden/compound-catalog-full.json` evaluates ten pure-property correlation families for all 14 records through DWSIM 9.0.5.0. Its normalized repeat capture is required to match before replacement.
- `tests/golden/compound-catalog-extended-equations.json` covers all 17 promoted records that use the additional constant, polynomial, exponential, or reduced-temperature equation forms; its normalized repeat also matches.
- Every JSON capture intended for parity is run twice and compared after removing `source.captured_utc`. `scripts/validate.py --compare FIRST SECOND` performs that normalized digest comparison and rejects property-read errors.
- The Phase 17 tank/PID tables are byte-identical executable exports. Phase 18 freezes official DWSIM clean-energy samples; the solar case was solved and saved again in DWSIM 10 before it was committed.

## Change control

1. Add or alter data only with an immutable source path or publication, source revision, units, and applicable validity range.
2. Recalculate equation vectors from the stated equation before hardcoding a test expectation.
3. Add loader rejection tests for malformed or semantically inconsistent records.
4. Capture DWSIM twice when executable parity is claimed and document any residual or model difference.
5. Version a schema when a change is not backward compatible; never reinterpret an existing field in place.

The three ChemSep correlation artifacts are generated together by
`scripts/extract_chemsep_correlations.py`. A no-argument run is the deterministic drift
gate; `--write` is the explicit regeneration operation.

The vendored DWSIM source is GPLv3. Public distribution of derived code or data requires a separate license review; the current 0.1.0 source promotion is for internal use and is not a public release tag.
