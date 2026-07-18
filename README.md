# me-sim-wrap

A process-simulation workspace with two parts:

1. **`src/mesim/`** — a cross-platform Python calculation kernel and HTTP API. This is the
   active, forward-looking implementation: it runs natively on macOS arm64 and Linux
   amd64/arm64 and never depends on the DWSIM runtime.
2. **`dwsim-windows/`** — the original Windows-only HTTP wrapper (`DWSIM.Api`) around the
   vendored DWSIM engine, kept as a reference oracle and legacy backend. See
   [`dwsim-windows/DWSIM.Api/README.md`](dwsim-windows/DWSIM.Api/README.md).

The Python kernel is the primary product. DWSIM is used only to generate deterministic
golden reference cases on Windows; it is **not** a runtime dependency of the Python service.

## Repository layout

```text
me-sim-wrap/
├── src/mesim/            ← Python kernel (the active implementation)
│   ├── api.py            ← FastAPI HTTP boundary
│   ├── compounds.py      ← versioned compound records
│   ├── streams.py        ← immutable stream / phase state + flash
│   ├── units.py          ← explicit unit-conversion boundary (SI internally)
│   ├── flowsheet.py      ← deterministic acyclic flowsheet execution
│   ├── thermo/           ← ideal, peng_robinson, flash, ...
│   └── unitops/          ← basic (mixer/valve/heater/...), pressure, ...
├── data/                 ← versioned compounds, interactions, correlations (JSON)
├── docs/                 ← compatibility matrix, plans, deployment notes
├── scripts/
│   ├── capture_dwsim_reference.ps1   ← Windows-only DWSIM golden-case capture
│   └── validate.py                   ← golden-case validation
├── tests/                ← unittest suites + tests/golden/ reference cases
├── dwsim-windows/        ← vendored DWSIM source + DWSIM.Api backend (Windows)
├── Dockerfile            ← Linux Python service image (no DWSIM)
└── pyproject.toml
```

## The Python kernel (`mesim`)

A deterministic process-simulation kernel with steady-state and bounded dynamic
capabilities. The original verified domain is five light hydrocarbons (methane,
ethane, propane, n-butane, n-pentane) with Peng-Robinson thermodynamics; later
golden-backed slices add the compounds required by reactions and columns.

**Design rules (from the rewrite plan):**

- Every input and output carries an explicit unit. Submitted value/unit are preserved;
  an SI value is derived for calculation and never overwrites the source.
- Compound constants, correlations, and binary interactions are immutable, versioned, and
  source-backed with provenance.
- Calculations are pure and deterministic. Solver failures are structured results
  (`converged`, `iterations`, `residual`, `algorithm`, `warnings`), not crashes.
- The kernel never imports or loads a DWSIM assembly.

### What is implemented

| Area | Status |
|---|---|
| Units boundary (temperature, pressure, flow, energy, power, enthalpy, density, transport, …) | Done |
| Five-compound ideal properties (cp, h, s, density, vapor pressure) | Done |
| Peng-Robinson pure + mixture (cubic roots, fugacity, departure h/s, stable-root selection) | Done (T1) |
| TP / PH / bubble / dew flash | Done (T2; documented DWSIM solver/model differences) |
| Activity-coefficient liquid VLE | Partial (T4: saved-source acetone/methanol NRTL bubble/dew slice) |
| Streams, mixer, splitter, heater, cooler, valve, equilibrium separator | Done (U0) |
| Acyclic flowsheet execution | Done |
| HTTP API (`/v1/*`) | Done (internal alpha) |
| Hydraulics | Partial (U3: pipes, defined heat-load and constant/gradient/tabulated defined-HTC plus estimated air/water/soil-HTC liquid thermal profiles with PR enthalpy coupling, optional insulation and local/global irradiation, fittings, orifice, two-phase correlations, relief sizing) |
| Reactors | Partial (U4: conversion, simultaneous vapor equilibrium reactions, and vapor Gibbs minimization) |
| Columns | Partial (U6/U7: shortcut, fixed-K absorber, fixed-profile material/energy gates, and live NRTL stage K/bubble-point parity) |
| Dynamics and controls | Partial (U9: holdup balances, fixed/adaptive ODE paths, tank/PID DWSIM parity, lumped HX, dynamic CSTR) |
| Specialty energy | Partial (U10: solar panel, wind turbine, and hydroelectric turbine source-equation parity) |

See [`docs/compatibility.md`](docs/compatibility.md) for the full matrix and current
parity/difference notes. Release boundaries are collected in
[`docs/model-limitations.md`](docs/model-limitations.md), with source lineage in
[`docs/data-provenance.md`](docs/data-provenance.md) and changes in
[`CHANGELOG.md`](CHANGELOG.md).

### HTTP API

The service is stateless and exposed via FastAPI. Requests carry explicit quantities
(`{ "value": ..., "unit": ... }`); responses echo inputs plus SI values.

| Method | Route | Purpose |
|---|---|---|
| GET | `/health` | Liveness + `schema_version` |
| GET | `/v1/compounds/{id}` | Source-backed compound record |
| POST | `/v1/flash/tp` | TP flash (phase, vapor fraction, compositions) |
| POST | `/v1/unitops/heater` | Heater to specified outlet temperature → outlet + duty |
| POST | `/v1/unitops/valve` | Isenthalpic valve to specified outlet pressure |
| POST | `/v1/flowsheets/u0` | Mixer → heater → valve → separator flowsheet |

```bash
curl -s -XPOST localhost:8000/v1/flash/tp \
  -H 'Content-Type: application/json' \
  -d '{"compound_ids":["C1","C2"],"composition":[0.5,0.5],
       "temperature":{"value":300,"unit":"K"},"pressure":{"value":10,"unit":"bar"}}'
```

Each calculation runs in a spawned worker with a timeout (`CALCULATION_TIMEOUT_S`, default
5 s) and a request-body cap (`MAX_REQUEST_BYTES`, default 1 MiB).

## Build, run, and test (Python)

Requires Python 3.12.

```bash
# install
pip install -e .

# run the unit/integration suites
python -m unittest discover -s tests -v

# validate golden reference cases (exit 0 only when all pass)
python scripts/validate.py --quiet

# run the API locally
uvicorn mesim.api:app --host 0.0.0.0 --port 8000
# or: fastapi dev  (see pyproject for entry point)
```

### Container

```bash
docker build -t me-sim .
docker run --rm -p 8000:8000 me-sim
```

The image is Linux-only, multi-arch (amd64/arm64), and **excludes `dwsim-windows/`**.

## Reference capture (Windows only)

Golden cases are produced from the vendored DWSIM engine on x64 Windows via
`scripts/capture_dwsim_reference.ps1`, then compared for determinism with
`scripts/validate.py --compare`. See `docs/compatibility.md` and
`docs/plans/2026-07-13-python-process-simulator-rewrite.md`.

## Deployment & security

- The service is **private-network only**: no application auth, must not be exposed on a
  public IP/ingress. Add bearer-token auth before any broader exposure. See
  [`docs/internal-alpha.md`](docs/internal-alpha.md).
- Required commands and acceptance tolerances are defined in the rewrite plan
  (`docs/plans/2026-07-13-python-process-simulator-rewrite.md`).

## License note

DWSIM is GPLv3 — upstream: [github.com/DanWBS/dwsim](https://github.com/DanWBS/dwsim),
license: [GNU GPL v3](https://www.gnu.org/licenses/gpl-3.0.html).
`dwsim-windows/` is vendored for reference capture only. Any future distribution of
DWSIM-derived code must be reviewed for GPLv3 compliance. The `mesim` Python package is
independent of the DWSIM runtime.
