# me-sim-wrap

Web-wrapping app around the **DWSIM** chemical process simulator. Goal: expose the DWSIM
engine over HTTP and build a browser frontend on top of it — **without rewriting the engine**.

## Layout

```
me-sim-wrap/
├─ CLAUDE.md              ← this file
├─ PLAN.md                ← phased plan + status (read this on restart)
└─ dwsim-windows/         ← DWSIM source (vendored upstream) + our backend
   ├─ DWSIM.Api/          ← OUR CODE: HTTP wrapper over the engine (Automation3)
   ├─ DWSIM.Automation/   ← engine entry point we call
   ├─ DWSIM.Thermodynamics, DWSIM.UnitOperations, DWSIM.FlowsheetSolver, ... (engine)
   └─ DWSIM.UI.Desktop.*  ← desktop UIs — we do NOT touch these
```

`DWSIM.Api/` lives _inside_ `dwsim-windows/` on purpose: its `.csproj` uses relative
`..\DWSIM.Automation` project references. Keeping it there means the existing `DWSIM.sln`
and relative refs resolve with zero rewiring. It is our code, not upstream.

## The one constraint that shapes everything

The engine targets **.NET Framework 4.6.2 → Windows only**, and loads **native x64 DLLs**
(CoolProp, Reaktoro). Consequences:

- Build + run only on **Windows** with MSBuild / Visual Studio. Not macOS, not Linux.
- Docker = **Windows containers on a Windows host** only. Won't run on a Mac. See `DWSIM.Api/README.md`.
- Do **not** port the engine to JS or .NET Core for a prototype — years of work + numeric risk.
  Keep it as a Windows backend service; put the web UI in front of it.

## Architecture

```
browser (web UI)  ──HTTP──>  DWSIM.Api (self-hosted Web API, .exe)  ──in-proc──>  DWSIM engine
```

The whole engine API surface we need is ~6 calls on `Automation3`:
`CreateFlowsheet` / `LoadFlowsheet2` / `CalculateFlowsheet4` / `SaveFlowsheet2`, plus
`IFlowsheet.GetFlowsheetSimulationObject(tag)` and `ISimulationObject.Get/SetPropertyValue`.

## Build & run (Windows)

```
cd dwsim-windows
nuget restore DWSIM.sln
msbuild DWSIM.Api\DWSIM.Api.csproj /p:Configuration=Release /p:Platform=x64
DWSIM.Api\bin\x64\Release\DWSIM.Api.exe        # serves http://localhost:9000
```

If VS/sln doesn't list the project: `dotnet sln DWSIM.sln add DWSIM.Api\DWSIM.Api.csproj`.
Platform (x64/x86) must match the engine's native deps — default is x64.

## Rules for working here

- Engine code (`dwsim-windows/DWSIM.*` except `DWSIM.Api`) is **vendored upstream** — avoid editing;
  if you must, note it, because it complicates pulling upstream updates.
- All new backend code goes in `DWSIM.Api/`. Frontend goes in a new top-level `web/` (Phase 2).
- The engine is **not thread-safe per flowsheet** — serialize solve calls per session (already done
  via a per-session lock in `FlowsheetStore`).
- No auth/CORS yet — add before exposing beyond localhost.

## Status

Phase 0 (backend scaffold) done — code is written but **never built/run** (was authored on macOS).
First Windows task: restore + build + run, confirm the engine boots and one flowsheet solves.
See `PLAN.md`.
