import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim.compounds import load_compounds, load_pr_interactions
from mesim.errors import OutOfRangeError, ValidationError
from mesim.thermo.ideal import load_correlations
from mesim.thermo.pure import load_saturated_liquid_correlations
from mesim.thermo.systems import PengRobinsonSystem
from mesim.thermo.transport import load_transport_correlations


ROOT = Path(__file__).parents[1]


class PureComponentCorrelationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.compounds = {
            record.id: record
            for record in load_compounds(ROOT / "data/compounds/v1.json")
        }
        cls.ideal = {
            record.compound_id: record
            for record in load_correlations(ROOT / "data/correlations/ideal-v1.json")
        }
        cls.transport = {
            record.compound_id: record
            for record in load_transport_correlations(
                ROOT / "data/correlations/transport-v1.json"
            )
        }
        cls.saturated = {
            record.compound_id: record
            for record in load_saturated_liquid_correlations(
                ROOT / "data/correlations/saturated-liquid-v1.json"
            )
        }

    def test_all_catalog_compounds_have_all_extracted_property_families(self):
        expected = set(self.compounds)
        self.assertEqual(len(expected), 391)
        self.assertEqual(set(self.ideal), expected)
        self.assertEqual(set(self.transport), expected)
        self.assertEqual(set(self.saturated), expected)
        result = subprocess.run(
            [sys.executable, "scripts/extract_chemsep_correlations.py"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_every_extracted_correlation_has_a_positive_midpoint_value(self):
        for record in self.ideal.values():
            for correlation, evaluate in (
                (record.heat_capacity_correlation, record.heat_capacity),
                (record.vapor_pressure_correlation, record.vapor_pressure),
            ):
                evaluate((correlation.minimum_k + correlation.maximum_k) / 2.0)
        for record in self.saturated.values():
            for correlation, evaluate in (
                (record.liquid_density_correlation, record.liquid_molar_density),
                (record.liquid_heat_capacity_correlation, record.liquid_heat_capacity),
                (record.heat_of_vaporization_correlation, record.heat_of_vaporization),
                (record.surface_tension_correlation, record.surface_tension),
            ):
                evaluate((correlation.minimum_k + correlation.maximum_k) / 2.0)
        for record in self.transport.values():
            for correlation in (
                record.liquid_viscosity,
                record.vapor_viscosity,
                record.liquid_thermal_conductivity,
                record.vapor_thermal_conductivity,
            ):
                correlation.value(
                    (correlation.minimum_k + correlation.maximum_k) / 2.0
                )

    def test_extracted_properties_match_dwsim_catalog_evaluation(self):
        captured = json.loads(
            (ROOT / "tests/golden/compound-catalog-full.json").read_text(
                encoding="utf-8-sig"
            )
        )["inputs"]["compounds"]
        self.assertEqual(len(captured), 14)
        for record in captured:
            compound_id = record["id"]
            molecular_weight = self.compounds[compound_id].molecular_weight.value
            ideal = self.ideal[compound_id]
            pure = self.saturated[compound_id]
            transport = self.transport[compound_id]
            ideal_reference = record["ideal_reference"]
            pure_reference = record["pure_reference"]

            heat_capacity = ideal.heat_capacity(
                ideal_reference["heat_capacity_temperature"]["value"], True
            ).value / molecular_weight / 1000.0
            vapor_pressure = ideal.vapor_pressure(
                ideal_reference["vapor_pressure_temperature"]["value"], True
            ).value
            self.assertTrue(
                math.isclose(
                    heat_capacity,
                    ideal_reference["heat_capacity"]["value"],
                    rel_tol=2.0e-15,
                )
            )
            self.assertTrue(
                math.isclose(
                    vapor_pressure,
                    ideal_reference["vapor_pressure"]["value"],
                    rel_tol=2.0e-15,
                )
            )

            values = {
                "liquid_density": pure.liquid_molar_density(
                    pure_reference["liquid_density"]["temperature"]["value"]
                ).value
                * molecular_weight,
                "liquid_heat_capacity": pure.liquid_heat_capacity(
                    pure_reference["liquid_heat_capacity"]["temperature"]["value"]
                ).value
                / molecular_weight
                / 1000.0,
                "heat_of_vaporization": pure.heat_of_vaporization(
                    pure_reference["heat_of_vaporization"]["temperature"]["value"]
                ).value
                / molecular_weight
                / 1000.0,
                "surface_tension": pure.surface_tension(
                    pure_reference["surface_tension"]["temperature"]["value"]
                ).value,
                "liquid_viscosity": transport.liquid_viscosity.value(
                    pure_reference["liquid_viscosity"]["temperature"]["value"], True
                ),
                "vapor_viscosity": transport.vapor_viscosity.value(
                    pure_reference["vapor_viscosity"]["temperature"]["value"], True
                ),
                "liquid_thermal_conductivity": transport.liquid_thermal_conductivity.value(
                    pure_reference["liquid_thermal_conductivity"]["temperature"]["value"],
                    True,
                ),
                "vapor_thermal_conductivity": transport.vapor_thermal_conductivity.value(
                    pure_reference["vapor_thermal_conductivity"]["temperature"]["value"],
                    True,
                ),
            }
            for name, value in values.items():
                self.assertTrue(
                    math.isclose(
                        value,
                        pure_reference[name]["value"]["value"],
                        rel_tol=2.0e-15,
                    ),
                    f"{compound_id} {name}",
                )

    def test_saturated_properties_enforce_range_and_equation_domain(self):
        methane = self.saturated["Methane"]
        with self.assertRaises(OutOfRangeError):
            methane.surface_tension(300.0)
        with self.assertRaises(ValidationError):
            methane.heat_of_vaporization(
                methane.critical_temperature_k, allow_extrapolation=True
            )
        with self.assertRaises(ValidationError):
            methane.liquid_molar_density(float("nan"))

    def test_loader_rejects_unknown_schema_and_invalid_equation(self):
        source = json.loads(
            (ROOT / "data/correlations/saturated-liquid-v1.json").read_text()
        )
        for mutation in ("schema", "equation"):
            broken = json.loads(json.dumps(source))
            if mutation == "schema":
                broken["schema_version"] = "unknown"
            else:
                broken["correlations"][0]["liquid_density"]["equation"] = 101
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "pure.json"
                path.write_text(json.dumps(broken))
                with self.assertRaises(ValidationError):
                    load_saturated_liquid_correlations(path)

    def test_pr_system_exposes_full_pure_component_records(self):
        methane = self.compounds["Methane"]
        system = PengRobinsonSystem(
            (methane,),
            load_pr_interactions(ROOT / "data/interactions/pr-v1.json"),
            tuple(self.ideal.values()),
            tuple(self.transport.values()),
            tuple(self.saturated.values()),
        )
        self.assertIs(system.ideal("Methane"), self.ideal["Methane"])
        self.assertIs(system.transport("Methane"), self.transport["Methane"])
        self.assertIs(
            system.saturated_liquid("Methane"), self.saturated["Methane"]
        )
        with self.assertRaises(ValidationError):
            system.saturated_liquid("Acetone")


if __name__ == "__main__":
    unittest.main()
