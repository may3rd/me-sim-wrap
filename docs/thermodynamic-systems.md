# Thermodynamic systems

`mesim.thermo.systems` binds ordered compounds and versioned model data into explicit,
immutable calculation systems. Unit operations consume a system capability instead of
loading correlation files or selecting model equations themselves.

| Stable model ID | Class | Verified capabilities | Explicit boundary |
|---|---|---|---|
| `peng-robinson-classic` | `PengRobinsonSystem` | TP, PH, and PS flash, total flash enthalpy, and source-backed ideal, saturated-liquid, and transport records | Classic PR with supplied binary interactions; a mixture is accepted only when every requested pair is explicit |
| `ideal-raoult` | `IdealRaoultSystem` | Liquid/vapor fugacity coefficients, K-values, bubble/dew pressure, and TP flash | DWSIM Ideal/Raoult equilibrium contract over an ordered source-backed vapor-pressure domain |
| `soave-redlich-kwong` | `SoaveRedlichKwongSystem` | Liquid/vapor/stable phase states and VLE TP flash | Classic SRK with its own frozen DWSIM binary-interaction table; PR interaction records are rejected |
| `peng-robinson-1978` | `PengRobinson1978System` | Liquid/vapor/stable phase states | PR78's piecewise alpha law, including the high-acentric-factor branch, with the frozen PR interaction domain |
| `peng-robinson-lee-kesler` | `PengRobinsonLeeKeslerSystem` | Liquid/vapor/stable phase states and VLE TP flash | DWSIM PR/LK's classic-PR equilibrium path; Lee-Kesler caloric and compressibility-property overrides remain explicitly outside this boundary |
| `peng-robinson-stryjek-vera-2-margules` | `PRSV2MargulesSystem` | Liquid/vapor/stable phase states | DWSIM PRSV2 alpha correction plus its composition-dependent asymmetric Margules mixing rule; 90 alpha records and eight interaction pairs are frozen from the source assets |
| `peng-robinson-stryjek-vera-2-van-laar` | `PRSV2VanLaarSystem` | Liquid/vapor/stable phase states | DWSIM PRSV2 alpha correction plus its rational asymmetric Van Laar mixing rule and separate eight-pair source table |
| `peng-robinson-1978-advanced` | `PengRobinson1978AdvancedSystem` | Liquid/vapor/stable phase states | PR78 plus all 13 DWSIM advanced mercury interaction expressions evaluated as source-normalized temperature polynomials; absent/zero expressions fall back to ordinary PR interactions |
| `soave-redlich-kwong-advanced` | `SoaveRedlichKwongAdvancedSystem` | Liquid/vapor/stable phase states and VLE TP flash | SRK plus a saved-case T/P interaction-expression dictionary; the parity case configures the DWSIM Mercury/N-pentane temperature polynomial before capture |
| `nrtl-acetone-methanol` | `NRTLSystem` | Modified-Raoult equilibrium ratios, bubble/dew pressure, bubble temperature, and Excess-mode phase enthalpies | Saved-source binary acetone/methanol domain only; no general TP/PH flash |
| `wilson-acetone-methanol` | `WilsonSystem` | Liquid activity coefficients | Complete DWSIM 9.0.5 Wilson interaction table with an exact 298.15 K molar-volume basis frozen only for acetone/methanol; no general flash or caloric model |
| `uniquac-1-propanol-water` | `UniquacSystem` | Liquid activity coefficients | Complete 376-pair DWSIM source table with resolved `R`, `Q`, and directional parameters only for 1-propanol/water; no general flash, LLE, or caloric model |
| `unifac-1-propanol-water` | `UnifacSystem` | Original-UNIFAC liquid activity coefficients | Complete installed 119-subgroup and 1,403-directed-interaction domains with a resolved 1-propanol/water group basis; no general flash, LLE, or caloric model |
| `unifac-ll-1-propanol-water` | `UnifacLLSystem` | UNIFAC-LL liquid activity coefficients | Shared 119-subgroup domain plus all 1,467 LLE-directed interactions and DWSIM's executable mixed-group normalization for 1-propanol/water; no general flash or LLE split solver |
| `modfac-dortmund-1-propanol-water` | `ModfacDortmundSystem` | Modified UNIFAC (Dortmund) liquid activity coefficients | Complete 108-subgroup and 1,167 paired temperature-dependent interaction domain with scoped 1-propanol/water group vectors; no general flash or caloric model |
| `modfac-nist-1-propanol-water` | `ModfacNistSystem` | Modified UNIFAC (NIST) liquid activity coefficients | Complete 201-subgroup and 1,969-directed temperature-dependent interaction domain with scoped 1-propanol/water group vectors; no general flash or caloric model |
| `chao-seader-methane-n-pentane` | `ChaoSeaderSystem` | Chao-Seader liquid/vapor fugacity coefficients and TP flash | Exact methane/N-pentane pure constants and source-equation parity at 350 K and 1 MPa; no Lee-Kesler calorics or broader compound domain |
| `grayson-streed-methane-n-pentane` | `GraysonStreedSystem` | Grayson-Streed liquid/vapor fugacity coefficients and TP flash | Exact methane/N-pentane pure constants and source-equation parity at 350 K and 1 MPa; no Lee-Kesler calorics or broader compound domain |

The registry is a fixed dictionary from stable model ID to constructor. Runtime plugin
registration and silent fallback are unsupported. An unknown ID, incomplete correlation
set, missing binary interaction, or compound-order mismatch fails before calculation.

PRSV2 data is regenerated by `scripts/extract_dwsim_prsv2_data.ps1` from the three
embedded DWSIM source assets, with each source hash retained in the generated JSON.
Repeat captures cover methane/ethane alpha behavior and the nonzero
acetone/cyclohexane Margules and Van Laar pairs. Feed-phase fugacity coefficients
agree within `5e-11` relative; TP flash and caloric methods remain outside these
first PRSV2 gates.

The advanced PR78/SRK source domain is separate from the 408-compound pure-property
domain: DWSIM's built-in override table contains only Mercury pairs. The scoped
`advanced-eos-v1.json` compound file freezes Mercury's EOS constants from `chemsep2.xml`,
while `prsrk-advanced-v1.json` preserves all 13 raw expressions and normalized
polynomial coefficients. This does not imply that Mercury has the full ideal,
transport, or saturated-liquid correlation capabilities of the main catalog.

Wilson data is regenerated from the installed DWSIM 9.0.5 embedded resource and a saved
acetone/methanol case by `scripts/extract_dwsim_wilson_data.ps1`. All 364 interaction
pairs are retained in cal/mol, while the two scoped molar volumes reproduce the package's
298.15 K `AUX_LIQDENSi` basis. The equation-level activity coefficients match both
repeat captures within floating-point precision. Fugacity, flash, and caloric behavior
remain outside this first Wilson system boundary.

UNIQUAC data is regenerated by `scripts/extract_dwsim_uniquac_data.ps1`. The artifact
retains every record from the embedded `uniquac.dat` resource and its hash, including
alternative regressions whose first source row wins in DWSIM, then freezes
the exact case-distinct names, ChemSep IDs, `R`/`Q` values, and resolved cal/mol parameters
for the 1-propanol/water parity domain. Repeatable DWSIM activity coefficients agree
within `2e-13` relative.

Original UNIFAC data is regenerated by `scripts/extract_dwsim_unifac_data.ps1` from
the installed runtime dictionaries. The source hashes, all subgroup `R`/`Q` records,
all directed primary-group energies, and the scoped compound surface-fraction vectors
are retained. The Python group-contribution equation matches repeat DWSIM activity
coefficients within floating-point precision.

UNIFAC-LL uses the same group-contribution equation with the separately hashed
`unifac_ll_ip.txt` matrix. Its DWSIM package builds the denominator from
`MODFACGroups` while `RET_VN` reads `UNIFACGroups`; the resulting non-unit surface
fraction sums are captured and preserved as executable behavior. The scoped activity
coefficients match the repeatable DWSIM vector exactly.

The Dortmund extractor freezes both embedded resource hashes, all `R`/`Q` groups,
and all paired six-coefficient interaction records. Its `r^(3/4)` combinatorial rule
and quadratic temperature interaction rule reproduce the saved activity vector within
`3e-15` relative.

The same extractor selects Modified UNIFAC (NIST) by exact case-distinct package name
and freezes its separate group and interaction resources. NIST stores its 1,969
directions as individual three-coefficient rows; direct-direction lookup reproduces
both repeat captures exactly for the scoped 1-propanol/water state.

Chao-Seader data is regenerated by `scripts/extract_dwsim_chao_seader_data.ps1`.
The artifact freezes the two case compounds' critical, acentric, liquid-volume, and
solubility constants together with runtime and vendored equation-source hashes. The
liquid regular-solution/pure-reference equation and vapor Redlich-Kwong equation match
repeat phase-fugacity captures within floating-point precision. The Python TP flash
converges the same equations more tightly than DWSIM's four-iteration reference; phase
fractions and compositions use a documented `1e-4` absolute reference tolerance.
The parameterized extractor also freezes Grayson-Streed's separate package/model source
hashes. Grayson-Streed shares the vapor Redlich-Kwong equation and compound constants,
but its distinct pure-liquid reference coefficients are evaluated and parity-tested
independently; its saved four-iteration flash uses the same reference tolerance.

The frozen compound catalog and all three pure-property datasets cover the same 408
case-distinct DWSIM/ChemSep names. These are the records among the 431-source catalog
whose ten required property families use equation types implemented by this runtime;
the other 23 remain explicitly excluded. The PR system exposes supported records through `ideal`,
`saturated_liquid`, and `transport`. This is property coverage, not a claim that all
possible binary combinations are supported: `pr-v1.json` contains 197 explicit pairs
(181 first-source records plus 16 accepted explicit-zero closures),
and the constructor rejects an unrecorded pair instead of assuming a zero interaction.
`chemsep-exclusions-v1.json` records the missing fields and unsupported equation IDs for
every excluded source compound.

The systems deliberately do not implement one artificial common flash interface. PR and
the accepted NRTL slice have different verified capabilities, and callers must request a
method the selected system actually provides. Adding another model requires its own
source-backed data domain and parity gate before its constructor is added to the registry.

The HTTP TP-flash boundary constructs `PengRobinsonSystem`. Live NRTL column solvers
consume `NRTLSystem`, which keeps activity and caloric data access out of the unit-operation
layer. The older three-data-argument `flash_stream` form remains temporarily accepted for
internal migration; new callers should pass a `PengRobinsonSystem` directly.

`scripts/extract_chemsep_correlations.py` deterministically regenerates the ideal,
transport, and saturated-liquid datasets from the vendored ChemSep XML. Run it without
arguments as a drift check or with `--write` after deliberately changing the frozen
source revision.
