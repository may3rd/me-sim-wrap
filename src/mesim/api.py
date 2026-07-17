import atexit
from dataclasses import asdict
from multiprocessing import TimeoutError as MultiprocessingTimeoutError, get_context
from pathlib import Path
import sys
from threading import Lock

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel, ConfigDict

from . import __version__
from .compounds import load_compounds, load_pr_interactions
from .errors import ValidationError
from .streams import PhaseState, StreamState, flash_stream
from .thermo.flash import tp_flash
from .thermo.ideal import load_correlations
from .units import Quantity
from .unitops.basic import equilibrium_separator, heater, mix_streams, valve


_SOURCE_DATA = Path(__file__).resolve().parents[2] / "data"
_INSTALLED_DATA = Path(sys.prefix) / "share" / "mesim" / "data"
DATA = _SOURCE_DATA if _SOURCE_DATA.is_dir() else _INSTALLED_DATA
COMPOUNDS = {compound.id: compound for compound in load_compounds(DATA / "compounds/v1.json")}
INTERACTIONS = load_pr_interactions(DATA / "interactions/pr-v1.json")
CORRELATIONS = load_correlations(DATA / "correlations/ideal-v1.json")
CALCULATION_TIMEOUT_S = 5.0
MAX_REQUEST_BYTES = 1_048_576
CALCULATION_WORKERS = 2

_CALCULATION_POOL = None
_CALCULATION_POOL_LOCK = Lock()

app = FastAPI(title="me-sim", version=__version__)


class _RequestTooLarge(Exception):
    pass


class RequestSizeLimit:
    def __init__(self, app, maximum: int) -> None:
        self.app = app
        self.maximum = maximum

    async def __call__(self, scope, receive, send) -> None:
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        length = dict(scope["headers"]).get(b"content-length")
        if length is not None and int(length) > self.maximum:
            await send({"type": "http.response.start", "status": 413, "headers": []})
            await send({"type": "http.response.body", "body": b'{"detail":"request body too large"}'})
            return
        received = 0

        async def limited_receive():
            nonlocal received
            message = await receive()
            if message["type"] == "http.request":
                received += len(message.get("body", b""))
                if received > self.maximum:
                    raise _RequestTooLarge
            return message

        try:
            await self.app(scope, limited_receive, send)
        except _RequestTooLarge:
            await send({"type": "http.response.start", "status": 413, "headers": []})
            await send({"type": "http.response.body", "body": b'{"detail":"request body too large"}'})


app.add_middleware(RequestSizeLimit, maximum=MAX_REQUEST_BYTES)


class QuantityInput(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)

    value: float
    unit: str

    def quantity(self, dimension: str) -> Quantity:
        return Quantity.from_value(self.value, self.unit, dimension)


class TPFlashRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)

    compound_ids: tuple[str, ...]
    composition: tuple[float, ...]
    temperature: QuantityInput
    pressure: QuantityInput


class StreamInput(TPFlashRequest):
    molar_flow: QuantityInput


class HeaterRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    stream: StreamInput
    outlet_temperature: QuantityInput


class ValveRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    stream: StreamInput
    outlet_pressure: QuantityInput
    temperature_bracket: tuple[QuantityInput, QuantityInput]


class U0FlowsheetRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    feeds: tuple[StreamInput, StreamInput]
    mixer_outlet_pressure: QuantityInput
    mixer_temperature_bracket: tuple[QuantityInput, QuantityInput]
    heater_outlet_temperature: QuantityInput
    valve_outlet_pressure: QuantityInput
    valve_temperature_bracket: tuple[QuantityInput, QuantityInput]


@app.exception_handler(ValidationError)
async def validation_error(_: Request, error: ValidationError) -> JSONResponse:
    return JSONResponse(status_code=422, content={"detail": str(error)})


def _calculation_result(function, request) -> tuple:
    try:
        return "ok", function(request)
    except HTTPException as error:
        return "http", error.status_code, error.detail
    except ValidationError as error:
        return "validation", str(error)


def _calculation_pool():
    global _CALCULATION_POOL
    with _CALCULATION_POOL_LOCK:
        if _CALCULATION_POOL is None:
            _CALCULATION_POOL = get_context("spawn").Pool(processes=CALCULATION_WORKERS)
        return _CALCULATION_POOL


def _discard_calculation_pool(pool) -> None:
    global _CALCULATION_POOL
    discard = False
    with _CALCULATION_POOL_LOCK:
        if _CALCULATION_POOL is pool:
            _CALCULATION_POOL = None
            discard = True
    if discard:
        pool.terminate()
        pool.join()


def _close_calculation_pool() -> None:
    global _CALCULATION_POOL
    with _CALCULATION_POOL_LOCK:
        pool = _CALCULATION_POOL
        _CALCULATION_POOL = None
    if pool is not None:
        pool.close()
        pool.join()


atexit.register(_close_calculation_pool)


def _limited(function, request) -> dict[str, object]:
    pool = _calculation_pool()
    pending = pool.apply_async(_calculation_result, (function, request))
    try:
        result = pending.get(CALCULATION_TIMEOUT_S)
    except MultiprocessingTimeoutError:
        _discard_calculation_pool(pool)
        raise HTTPException(status_code=408, detail="calculation time limit exceeded")
    except (EOFError, OSError) as error:
        _discard_calculation_pool(pool)
        raise HTTPException(status_code=500, detail="calculation worker exited without a result") from error
    if result[0] == "ok":
        return result[1]
    if result[0] == "http":
        raise HTTPException(status_code=result[1], detail=result[2])
    raise ValidationError(result[1])


def _quantity(quantity: Quantity) -> dict[str, float | str]:
    return {"value": quantity.value, "unit": quantity.unit, "si_value": quantity.si_value}


def _compounds(compound_ids: tuple[str, ...]):
    try:
        return tuple(COMPOUNDS[compound_id] for compound_id in compound_ids)
    except KeyError as error:
        raise HTTPException(status_code=404, detail=f"unknown compound: {error.args[0]}") from error


def _phase(stream: StreamInput) -> PhaseState:
    temperature = stream.temperature.quantity("temperature")
    pressure = stream.pressure.quantity("pressure")
    flow = stream.molar_flow.quantity("molar_flow")
    return flash_stream(
        StreamState(temperature.si_value, pressure.si_value, flow.si_value, stream.compound_ids, stream.composition),
        _compounds(stream.compound_ids), INTERACTIONS, CORRELATIONS,
    )


def _stream_input(stream: StreamInput) -> dict[str, object]:
    return {
        "compound_ids": list(stream.compound_ids), "composition": list(stream.composition),
        "temperature": _quantity(stream.temperature.quantity("temperature")),
        "pressure": _quantity(stream.pressure.quantity("pressure")),
        "molar_flow": _quantity(stream.molar_flow.quantity("molar_flow")),
    }


def _phase_response(phase: PhaseState) -> dict[str, object]:
    return {
        "stream": {
            "temperature_k": phase.stream.temperature_k, "pressure_pa": phase.stream.pressure_pa,
            "molar_flow_kmol_s": phase.stream.molar_flow_kmol_s, "compound_ids": phase.stream.compound_ids,
            "composition": phase.stream.composition,
        },
        "phase": phase.flash.phase, "vapor_fraction": phase.flash.vapor_fraction,
        "enthalpy_j_per_kmol": phase.enthalpy_j_per_kmol,
    }


def _stream_response(stream: StreamState) -> dict[str, object]:
    return {
        "temperature_k": stream.temperature_k, "pressure_pa": stream.pressure_pa,
        "molar_flow_kmol_s": stream.molar_flow_kmol_s, "compound_ids": stream.compound_ids,
        "composition": stream.composition, "enthalpy_j_per_kmol": stream.enthalpy_j_per_kmol,
    }


@app.get("/health")
def health() -> dict[str, str]:
    return {"schema_version": "mesim-api-1", "status": "ok"}


def _bracket(items: tuple[QuantityInput, QuantityInput]) -> tuple[Quantity, Quantity]:
    return tuple(item.quantity("temperature") for item in items)


@app.get("/v1/compounds/{compound_id}")
def compound(compound_id: str) -> dict[str, object]:
    try:
        record = COMPOUNDS[compound_id]
    except KeyError as error:
        raise HTTPException(status_code=404, detail=f"unknown compound: {compound_id}") from error
    return {"schema_version": "mesim-api-1", "compound": asdict(record)}


@app.post("/v1/flash/tp")
def flash(request: TPFlashRequest) -> dict[str, object]:
    return _limited(_flash, request)


def _flash(request: TPFlashRequest) -> dict[str, object]:
    compounds = _compounds(request.compound_ids)
    temperature = request.temperature.quantity("temperature")
    pressure = request.pressure.quantity("pressure")
    result = tp_flash(compounds, request.composition, INTERACTIONS, temperature.si_value, pressure.si_value)
    return {
        "schema_version": "mesim-api-1",
        "inputs": {
            "compound_ids": list(request.compound_ids), "composition": list(request.composition),
            "temperature": _quantity(temperature), "pressure": _quantity(pressure),
        },
        "result": {
            "phase": result.phase, "vapor_fraction": result.vapor_fraction,
            "liquid_composition": result.liquid_composition, "vapor_composition": result.vapor_composition,
            "report": asdict(result.report),
        },
    }


@app.post("/v1/unitops/heater")
def heat(request: HeaterRequest) -> dict[str, object]:
    return _limited(_heat, request)


def _heat(request: HeaterRequest) -> dict[str, object]:
    inlet = _phase(request.stream)
    outlet_temperature = request.outlet_temperature.quantity("temperature")
    result = heater(inlet, _compounds(request.stream.compound_ids), INTERACTIONS, CORRELATIONS, outlet_temperature.si_value)
    return {
        "schema_version": "mesim-api-1",
        "inputs": {"stream": _stream_input(request.stream), "outlet_temperature": _quantity(outlet_temperature)},
        "outlet": _phase_response(result.outlet),
        "energy": {"duty_w": result.energy.duty_w},
    }


@app.post("/v1/unitops/valve")
def throttle(request: ValveRequest) -> dict[str, object]:
    return _limited(_throttle, request)


def _throttle(request: ValveRequest) -> dict[str, object]:
    inlet = _phase(request.stream)
    outlet_pressure = request.outlet_pressure.quantity("pressure")
    bracket = _bracket(request.temperature_bracket)
    result = valve(
        inlet, _compounds(request.stream.compound_ids), INTERACTIONS, CORRELATIONS, outlet_pressure.si_value,
        tuple(item.si_value for item in bracket),
    )
    return {
        "schema_version": "mesim-api-1",
        "inputs": {
            "stream": _stream_input(request.stream), "outlet_pressure": _quantity(outlet_pressure),
            "temperature_bracket": [_quantity(item) for item in bracket],
        },
        "outlet": _phase_response(result),
    }


@app.post("/v1/flowsheets/u0")
def solve_u0(request: U0FlowsheetRequest) -> dict[str, object]:
    return _limited(_solve_u0, request)


def _solve_u0(request: U0FlowsheetRequest) -> dict[str, object]:
    first, second = (_phase(feed) for feed in request.feeds)
    compounds = _compounds(first.stream.compound_ids)
    mixer_pressure = request.mixer_outlet_pressure.quantity("pressure")
    mixer_bracket = _bracket(request.mixer_temperature_bracket)
    mixed = mix_streams((first, second), compounds, INTERACTIONS, CORRELATIONS, mixer_pressure.si_value, tuple(item.si_value for item in mixer_bracket))
    heater_temperature = request.heater_outlet_temperature.quantity("temperature")
    heated = heater(mixed, compounds, INTERACTIONS, CORRELATIONS, heater_temperature.si_value).outlet
    valve_pressure = request.valve_outlet_pressure.quantity("pressure")
    valve_bracket = _bracket(request.valve_temperature_bracket)
    throttled = valve(heated, compounds, INTERACTIONS, CORRELATIONS, valve_pressure.si_value, tuple(item.si_value for item in valve_bracket))
    separated = equilibrium_separator(throttled, compounds, CORRELATIONS)
    return {
        "schema_version": "mesim-api-1",
        "inputs": {
            "feeds": [_stream_input(feed) for feed in request.feeds],
            "mixer_outlet_pressure": _quantity(mixer_pressure),
            "mixer_temperature_bracket": [_quantity(item) for item in mixer_bracket],
            "heater_outlet_temperature": _quantity(heater_temperature),
            "valve_outlet_pressure": _quantity(valve_pressure),
            "valve_temperature_bracket": [_quantity(item) for item in valve_bracket],
        },
        "streams": {
            "mixer": _phase_response(mixed), "heater": _phase_response(heated), "valve": _phase_response(throttled),
            "liquid": _stream_response(separated.liquid) if separated.liquid else None,
            "vapor": _stream_response(separated.vapor) if separated.vapor else None,
        },
    }
