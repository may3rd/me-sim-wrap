# DWSIM.Api

Thin HTTP wrapper over the DWSIM engine (`Automation3`). Self-hosted console app — no IIS.
Windows + .NET Framework 4.6.2 only (the engine is Windows-bound). This is the backend for a web frontend.

## Build & run (Windows)

```
nuget restore DWSIM.sln          # or: msbuild /t:restore
msbuild DWSIM.Api\DWSIM.Api.csproj /p:Configuration=Release /p:Platform=x64
bin\x64\Release\DWSIM.Api.exe    # listens on http://localhost:9000
```

Add to the solution first if VS doesn't see it: `dotnet sln DWSIM.sln add DWSIM.Api\DWSIM.Api.csproj`
(or add an existing project in Visual Studio). Set `DWSIM_API_URL` to change the bind address.

> Platform must match the engine's native deps (CoolProp). This project defaults to **x64**; if your
> DWSIM build is x86, switch both. First `Automation3()` construction is slow — it loads all property packages.

## Endpoints

| Method | Route | Body | Purpose |
|---|---|---|---|
| POST | `/api/flowsheet` | — | New empty flowsheet → `{id}` |
| POST | `/api/flowsheet/upload` | raw `.dwxmz` bytes | Load a saved model → `{id}` |
| GET | `/api/flowsheet/{id}/objects` | — | List all objects (tag, type, calculated) |
| GET | `/api/flowsheet/{id}/objects/{tag}` | — | All properties + values for one object |
| POST | `/api/flowsheet/{id}/property` | `{tag, property, value}` | Set one input |
| POST | `/api/flowsheet/{id}/solve` | — | Run solver → `{success, errors[]}` |
| GET | `/api/flowsheet/{id}/download` | — | Download solved `.dwxmz` |
| DELETE | `/api/flowsheet/{id}` | — | Drop the session |

## Example

```bash
ID=$(curl -s -XPOST localhost:9000/api/flowsheet/upload --data-binary @model.dwxmz | jq -r .id)
curl -s -XPOST localhost:9000/api/flowsheet/$ID/property \
  -H 'Content-Type: application/json' -d '{"tag":"Feed","property":"PROP_MS_0","value":350}'
curl -s -XPOST localhost:9000/api/flowsheet/$ID/solve
curl -s localhost:9000/api/flowsheet/$ID/objects/Product
```

## Known limits (by design)

- In-memory session store, single process. Add Redis/DB for scale-out or persistence.
- One `lock` per flowsheet — the engine is not thread-safe. Concurrent solves on the *same* id serialize.
- No auth/CORS yet — add before exposing beyond localhost. Frontend on another origin needs CORS.
- Property ids are DWSIM's internal strings (`PROP_MS_0` = temperature, etc.). `GET .../objects/{tag}`
  lists the valid ones for each object.
```
