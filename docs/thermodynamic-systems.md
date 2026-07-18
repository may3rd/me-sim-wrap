# Thermodynamic systems

`mesim.thermo.systems` binds ordered compounds and versioned model data into explicit,
immutable calculation systems. Unit operations consume a system capability instead of
loading correlation files or selecting model equations themselves.

| Stable model ID | Class | Verified capabilities | Explicit boundary |
|---|---|---|---|
| `peng-robinson-classic` | `PengRobinsonSystem` | TP, PH, and PS flash, total flash enthalpy, and source-backed ideal, saturated-liquid, and transport records | Classic PR with supplied binary interactions; a mixture is accepted only when every requested pair is explicit |
| `nrtl-acetone-methanol` | `NRTLSystem` | Modified-Raoult equilibrium ratios, bubble/dew pressure, bubble temperature, and Excess-mode phase enthalpies | Saved-source binary acetone/methanol domain only; no general TP/PH flash |

The registry is a fixed dictionary from stable model ID to constructor. Runtime plugin
registration and silent fallback are unsupported. An unknown ID, incomplete correlation
set, missing binary interaction, or compound-order mismatch fails before calculation.

The frozen compound catalog and all three pure-property datasets cover the same 391
case-distinct DWSIM/ChemSep names. These are the records among the 431-source catalog
whose ten required property families use equation types implemented by this runtime;
the other 40 remain explicitly excluded. The PR system exposes supported records through `ideal`,
`saturated_liquid`, and `transport`. This is property coverage, not a claim that all
possible binary combinations are supported: `pr-v1.json` contains 196 explicit pairs
(180 first-source records plus 16 accepted explicit-zero closures),
and the constructor rejects an unrecorded pair instead of assuming a zero interaction.

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
