# Plan — me-sim-wrap

> The Windows-wrapper plan below is retained as historical context. The active rewrite plan is
> [`docs/plans/2026-07-13-python-process-simulator-rewrite.md`](docs/plans/2026-07-13-python-process-simulator-rewrite.md).

Wrap the DWSIM engine in a web app. Backend = thin HTTP layer over `Automation3` (Windows).
Frontend = browser UI. Build the cheap, useful parts first; the interactive canvas is last.

Legend: ✅ done · 🔜 next · ⬜ later

---

## Phase 0 — Backend scaffold ✅ (code only, unbuilt)

`DWSIM.Api/` self-hosted Web API console app. 8 endpoints: create/upload flowsheet,
list objects, read/set properties, solve, download, delete. In-memory session store,
per-session lock. All engine calls verified against the interfaces — but **never compiled**
(authored on macOS; engine needs Windows).

---

## Phase 1 — Build & verify on Windows 🔜 (START HERE)

Goal: prove the engine boots in our process and one flowsheet solves end-to-end.

1. `nuget restore DWSIM.sln` (or `msbuild /t:restore`).
2. Add project to solution if needed: `dotnet sln DWSIM.sln add DWSIM.Api\DWSIM.Api.csproj`.
3. Build: `msbuild DWSIM.Api\DWSIM.Api.csproj /p:Configuration=Release /p:Platform=x64`.
   - **Watch:** native DLLs (CoolProp/Reaktoro) must land next to `DWSIM.Api.exe`. If missing at
     runtime, copy them explicitly from the desktop build output or mark copy-local.
   - **Watch:** Newtonsoft.Json version conflict — `app.config` has a binding redirect; adjust the
     version if the engine ships a different one.
   - **Watch:** x64 vs x86 must match the native deps.
4. Run `DWSIM.Api.exe`. Confirm "Engine ready" prints (property packages loaded).
5. Verify with a real model: take a sample `.dwxmz` from DWSIM's samples, `POST /upload`,
   `POST /solve`, `GET /objects/{tag}` — confirm results are non-empty and match the desktop app.

**Exit criteria:** one known flowsheet solves via HTTP with numbers matching desktop DWSIM.

---

## Phase 2 — Minimal web frontend ⬜

New top-level `web/` (plain HTML/JS or a small SPA — decide at start of phase).
Forms + tables only, **no flowsheet canvas yet**:

- Upload a `.dwxmz`, list its objects.
- Select an object → show its properties (from `GET /objects/{tag}`).
- Edit an input property → `POST /property` → `POST /solve` → refresh results.
- Download the solved file.

**Exit criteria:** a non-developer can load, tweak, solve, and read results in a browser.

Add **CORS** to `DWSIM.Api` here (frontend will be a different origin during dev).

---

## Phase 3 — Interactive flowsheet canvas ⬜

The one large piece. Rebuild the desktop SkiaSharp drawing surface in the browser (SVG/Canvas):
render objects + connections, drag to place, connect ports. Solver stays server-side — this is
UI only. Scope carefully; consider read-only render first, then editing.

---

## Phase 4 — Harden & deploy ⬜

- Auth (the API is currently open).
- Session lifecycle: eviction/TTL (in-memory store grows unbounded today).
- Persistence: swap in-memory store for DB/Redis if multi-instance or durability needed.
- Deploy: Windows container on a Windows host, or a Windows VM/service. CI on `windows-latest`.

---

## Decisions parked for later

- **Frontend stack** (vanilla vs React) — decide at Phase 2.
- **Scale-out / multi-user** — only if needed; single-process is fine for dev.
- **.NET Core migration** — likely never; only if Linux/container hosting becomes a hard requirement,
  and only after confirming the native deps can follow.
