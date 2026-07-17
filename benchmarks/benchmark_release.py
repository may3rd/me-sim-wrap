"""Repeatable Phase 19 latency gate for representative kernel and API paths."""
from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fastapi.testclient import TestClient

from mesim.api import app
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.streams import StreamState, flash_stream
from mesim.thermo.flash import tp_flash
from mesim.thermo.ideal import load_correlations
from mesim.unitops.basic import valve
from mesim.unitops.dynamics import (
    PIDConfig,
    TankLevelControlConfig,
    simulate_dwsim_tank_level_control,
)


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    samples: int
    median_ms: float
    p95_ms: float
    maximum_ms: float
    threshold_p95_ms: float
    passed: bool


def _measure(function: Callable[[], object], samples: int, threshold_ms: float) -> BenchmarkResult:
    function()
    durations = []
    for _ in range(samples):
        started = time.perf_counter_ns()
        function()
        durations.append((time.perf_counter_ns() - started) / 1_000_000.0)
    durations.sort()
    p95_index = max(0, math_ceil(0.95 * samples) - 1)
    p95 = durations[p95_index]
    return BenchmarkResult(
        samples,
        statistics.median(durations),
        p95,
        max(durations),
        threshold_ms,
        p95 <= threshold_ms,
    )


def math_ceil(value: float) -> int:
    integer = int(value)
    return integer if value == integer else integer + 1


def _fixtures():
    compounds = {record.id: record for record in load_compounds(ROOT / "data/compounds/v1.json")}
    selected = (compounds["Methane"], compounds["Ethane"])
    interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
    correlations = load_correlations(ROOT / "data/correlations/ideal-v1.json")
    inlet = flash_stream(
        StreamState(190.0, 600_000.0, 2.0, ("Methane", "Ethane"), (0.7, 0.3)),
        selected,
        interactions,
        correlations,
    )
    density = 997.060396254973
    tank = TankLevelControlConfig(
        duration_s=600.0,
        step_s=5.0,
        tank_volume_m3=2.0,
        tank_height_m=2.0,
        liquid_density_kg_m3=density,
        inlet_mass_flow_kg_s=10.0,
        initial_contents_volume_m3=10.0 / density,
        initial_outlet_opening_percent=33.8440722851719,
        inlet_opening_percent=50.0,
        tank_base_pressure_pa=101325.0,
        downstream_pressure_pa=109909.260185509,
        gravity_m_s2=9.8,
        valve_kv=400.0,
        pid=PIDConfig(
            119.455297118019,
            4.52783924040604,
            16.3382185733144,
            1.7,
            1.7,
            20.0,
            True,
            -1000.0,
            1000.0,
        ),
    )
    return selected, interactions, correlations, inlet, tank


TP_PAYLOAD = {
    "compound_ids": ["Methane", "Ethane"],
    "composition": [0.7, 0.3],
    "temperature": {"value": 180.0, "unit": "K"},
    "pressure": {"value": 500.0, "unit": "kPa"},
}
VALVE_PAYLOAD = {
    "stream": {
        **TP_PAYLOAD,
        "temperature": {"value": 190.0, "unit": "K"},
        "pressure": {"value": 600.0, "unit": "kPa"},
        "molar_flow": {"value": 2.0, "unit": "kmol/s"},
    },
    "outlet_pressure": {"value": 300.0, "unit": "kPa"},
    "temperature_bracket": [
        {"value": 140.0, "unit": "K"},
        {"value": 240.0, "unit": "K"},
    ],
}
FEED = {
    "compound_ids": ["Methane", "Ethane"],
    "temperature": {"value": 180.0, "unit": "K"},
    "pressure": {"value": 500.0, "unit": "kPa"},
    "molar_flow": {"value": 1.0, "unit": "kmol/s"},
}
U0_PAYLOAD = {
    "feeds": [
        {**FEED, "composition": [1.0, 0.0]},
        {**FEED, "composition": [0.0, 1.0]},
    ],
    "mixer_outlet_pressure": {"value": 500.0, "unit": "kPa"},
    "mixer_temperature_bracket": [
        {"value": 140.0, "unit": "K"},
        {"value": 240.0, "unit": "K"},
    ],
    "heater_outlet_temperature": {"value": 190.0, "unit": "K"},
    "valve_outlet_pressure": {"value": 300.0, "unit": "kPa"},
    "valve_temperature_bracket": [
        {"value": 140.0, "unit": "K"},
        {"value": 240.0, "unit": "K"},
    ],
}


def run(core_samples: int, api_samples: int) -> dict[str, object]:
    selected, interactions, correlations, inlet, tank = _fixtures()
    client = TestClient(app)

    def checked_post(path: str, payload: dict[str, object]) -> None:
        response = client.post(path, json=payload)
        if response.status_code != 200:
            raise RuntimeError(f"benchmark request failed: {path} returned {response.status_code}")

    measurements = {
        "core_tp_flash": _measure(
            lambda: tp_flash(selected, (0.7, 0.3), interactions, 180.0, 500_000.0),
            core_samples,
            25.0,
        ),
        "core_ph_valve": _measure(
            lambda: valve(inlet, selected, interactions, correlations, 300_000.0, (140.0, 240.0)),
            max(5, core_samples // 4),
            500.0,
        ),
        "core_tank_121_points": _measure(
            lambda: simulate_dwsim_tank_level_control(tank),
            core_samples,
            10.0,
        ),
        "api_tp_flash": _measure(
            lambda: checked_post("/v1/flash/tp", TP_PAYLOAD),
            api_samples,
            1000.0,
        ),
        "api_ph_valve": _measure(
            lambda: checked_post("/v1/unitops/valve", VALVE_PAYLOAD),
            api_samples,
            1500.0,
        ),
        "api_u0_flowsheet": _measure(
            lambda: checked_post("/v1/flowsheets/u0", U0_PAYLOAD),
            api_samples,
            2500.0,
        ),
    }
    return {
        "schema_version": "mesim-benchmark-1",
        "environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
        },
        "basis": "single-process sequential local TestClient; p95 thresholds are below 50% of the 5 s API deadline",
        "measurements": {name: asdict(result) for name, result in measurements.items()},
        "native_acceleration_required": not all(result.passed for result in measurements.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--core-samples", type=int, default=40)
    parser.add_argument("--api-samples", type=int, default=5)
    args = parser.parse_args()
    if args.core_samples < 5 or args.api_samples < 3:
        parser.error("use at least five core samples and three API samples")
    report = run(args.core_samples, args.api_samples)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 1 if report["native_acceleration_required"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
