from dataclasses import asdict
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel, ConfigDict

from .compounds import load_compounds, load_pr_interactions
from .errors import ValidationError
from .thermo.flash import tp_flash
from .units import Quantity


DATA = Path(__file__).resolve().parents[2] / "data"
COMPOUNDS = {compound.id: compound for compound in load_compounds(DATA / "compounds/v1.json")}
INTERACTIONS = load_pr_interactions(DATA / "interactions/pr-v1.json")

app = FastAPI(title="me-sim", version="0.1.0a0")


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


@app.exception_handler(ValidationError)
async def validation_error(_: Request, error: ValidationError) -> JSONResponse:
    return JSONResponse(status_code=422, content={"detail": str(error)})


def _quantity(quantity: Quantity) -> dict[str, float | str]:
    return {"value": quantity.value, "unit": quantity.unit, "si_value": quantity.si_value}


@app.get("/v1/compounds/{compound_id}")
def compound(compound_id: str) -> dict[str, object]:
    try:
        record = COMPOUNDS[compound_id]
    except KeyError as error:
        raise HTTPException(status_code=404, detail=f"unknown compound: {compound_id}") from error
    return {"schema_version": "mesim-api-1", "compound": asdict(record)}


@app.post("/v1/flash/tp")
def flash(request: TPFlashRequest) -> dict[str, object]:
    try:
        compounds = tuple(COMPOUNDS[compound_id] for compound_id in request.compound_ids)
    except KeyError as error:
        raise HTTPException(status_code=404, detail=f"unknown compound: {error.args[0]}") from error
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
