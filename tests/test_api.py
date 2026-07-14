import math
import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fastapi.testclient import TestClient

from mesim.api import app


class ApiTest(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(app)

    def test_compound_lookup_returns_versioned_data(self):
        response = self.client.get("/v1/compounds/Methane")

        self.assertEqual(response.status_code, 200)
        body = response.json()
        self.assertEqual((body["schema_version"], body["compound"]["id"]), ("mesim-api-1", "Methane"))

    def test_tp_flash_preserves_submitted_units_and_returns_phase(self):
        response = self.client.post("/v1/flash/tp", json={
            "compound_ids": ["Methane", "Ethane"],
            "composition": [0.7, 0.3],
            "temperature": {"value": -93.15, "unit": "degC"},
            "pressure": {"value": 500.0, "unit": "kPa"},
        })

        self.assertEqual(response.status_code, 200)
        body = response.json()
        self.assertEqual(body["inputs"]["temperature"]["value"], -93.15)
        self.assertEqual(body["inputs"]["temperature"]["unit"], "degC")
        self.assertTrue(math.isclose(body["inputs"]["temperature"]["si_value"], 180.0, abs_tol=1e-12))
        self.assertEqual(body["inputs"]["pressure"], {"value": 500.0, "unit": "kPa", "si_value": 500_000.0})
        self.assertEqual(body["result"]["phase"], "two-phase")
        self.assertTrue(body["result"]["report"]["converged"])

    def test_api_rejects_unknown_compounds_and_units(self):
        self.assertEqual(self.client.get("/v1/compounds/Unknown").status_code, 404)
        response = self.client.post("/v1/flash/tp", json={
            "compound_ids": ["Methane"], "composition": [1.0],
            "temperature": {"value": 300.0, "unit": "banana"},
            "pressure": {"value": 1.0, "unit": "bar"},
        })
        self.assertEqual(response.status_code, 422)

    def test_heater_endpoint_returns_calculated_signed_duty(self):
        response = self.client.post("/v1/unitops/heater", json={
            "stream": {
                "compound_ids": ["Methane", "Ethane"], "composition": [0.7, 0.3],
                "temperature": {"value": 180.0, "unit": "K"},
                "pressure": {"value": 500.0, "unit": "kPa"},
                "molar_flow": {"value": 2.0, "unit": "kmol/s"},
            },
            "outlet_temperature": {"value": 190.0, "unit": "K"},
        })

        self.assertEqual(response.status_code, 200)
        body = response.json()
        self.assertEqual(body["outlet"]["stream"]["temperature_k"], 190.0)
        self.assertGreater(body["energy"]["duty_w"], 0.0)

    def test_valve_endpoint_isenthalpically_reduces_pressure(self):
        response = self.client.post("/v1/unitops/valve", json={
            "stream": {
                "compound_ids": ["Methane", "Ethane"], "composition": [0.7, 0.3],
                "temperature": {"value": 190.0, "unit": "K"},
                "pressure": {"value": 600.0, "unit": "kPa"},
                "molar_flow": {"value": 2.0, "unit": "kmol/s"},
            },
            "outlet_pressure": {"value": 300.0, "unit": "kPa"},
            "temperature_bracket": [{"value": 140.0, "unit": "K"}, {"value": 240.0, "unit": "K"}],
        })

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["outlet"]["stream"]["pressure_pa"], 300_000.0)


if __name__ == "__main__":
    unittest.main()
