# Thermodynamic systems

`mesim.thermo.systems` binds ordered compounds and versioned model data into explicit,
immutable calculation systems. Unit operations consume a system capability instead of
loading correlation files or selecting model equations themselves.

| Stable model ID | Class | Verified capabilities | Explicit boundary |
|---|---|---|---|
| `peng-robinson-classic` | `PengRobinsonSystem` | TP, PH, and PS flash plus total flash enthalpy | Classic PR with supplied binary interactions and ideal correlations |
| `nrtl-acetone-methanol` | `NRTLSystem` | Modified-Raoult equilibrium ratios, bubble/dew pressure, bubble temperature, and Excess-mode phase enthalpies | Saved-source binary acetone/methanol domain only; no general TP/PH flash |

The registry is a fixed dictionary from stable model ID to constructor. Runtime plugin
registration and silent fallback are unsupported. An unknown ID, incomplete correlation
set, missing binary interaction, or compound-order mismatch fails before calculation.

The systems deliberately do not implement one artificial common flash interface. PR and
the accepted NRTL slice have different verified capabilities, and callers must request a
method the selected system actually provides. Adding another model requires its own
source-backed data domain and parity gate before its constructor is added to the registry.

The HTTP TP-flash boundary constructs `PengRobinsonSystem`. Live NRTL column solvers
consume `NRTLSystem`, which keeps activity and caloric data access out of the unit-operation
layer. The older three-data-argument `flash_stream` form remains temporarily accepted for
internal migration; new callers should pass a `PengRobinsonSystem` directly.
