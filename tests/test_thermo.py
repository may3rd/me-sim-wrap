import json
import math
import sys
import unittest
from xml.etree import ElementTree
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import OutOfRangeError, ValidationError
from mesim.thermo.ideal import ideal_gas_density, load_correlations


ROOT = Path(__file__).parents[1]


class IdealPropertiesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.correlations = {c.compound_id: c for c in load_correlations(ROOT / "data/correlations/ideal-v1.json")}

    def test_heat_capacity_and_vapor_pressure_match_chemsep_equations(self):
        expected_cp = {
            "Methane": 35896.79873540413,
            "Ethane": 53054.54046190328,
            "Propane": 75253.45426427145,
            "N-butane": 101725.78477833504,
            "N-pentane": 120847.8186718621,
        }
        for compound_id, expected in expected_cp.items():
            self.assertTrue(math.isclose(self.correlations[compound_id].heat_capacity(300).value, expected, rel_tol=1e-12))
        self.assertTrue(math.isclose(self.correlations["N-pentane"].vapor_pressure(300).value, 72902.6246331062, rel_tol=1e-12))

    def test_enthalpy_entropy_and_density_obey_ideal_gas_relations(self):
        methane = self.correlations["Methane"]
        self.assertEqual(methane.enthalpy_change(298.15, 298.15).value, 0.0)
        self.assertEqual(methane.entropy_change(298.15, 101325, 298.15, 101325).value, 0.0)
        dh = methane.enthalpy_change(300.01, 300).value
        self.assertTrue(math.isclose(dh / 0.01, methane.heat_capacity(300.005).value, rel_tol=1e-8))
        pressure_entropy = methane.entropy_change(300, 202650, 300, 101325).value
        self.assertTrue(math.isclose(pressure_entropy, -8314.46261815324 * math.log(2), rel_tol=1e-12))
        density = ideal_gas_density(16.04246, 300, 101325)
        self.assertEqual(density.unit, "kg/m3")
        self.assertTrue(math.isclose(density.value, 0.6516766162577914, rel_tol=1e-12))
        self.assertEqual(methane.enthalpy_change(300, 298.15).unit, "J/kmol")
        self.assertEqual(methane.entropy_change(300, 101325).unit, "J/kmol/K")
        full_range = methane.entropy_change(1500, 101325, 10, 101325).value
        self.assertTrue(math.isclose(full_range, 207430.066745, rel_tol=1e-9))
        for pressure in (math.nan, math.inf):
            with self.assertRaises(ValidationError):
                methane.entropy_change(300, pressure)
        extreme = methane.entropy_change(300, 5e-324, 300, 1e308).value
        self.assertTrue(math.isfinite(extreme))

    def test_ranges_are_enforced_and_extrapolation_is_flagged(self):
        methane = self.correlations["Methane"]
        with self.assertRaises(OutOfRangeError):
            methane.vapor_pressure(300)
        result = methane.vapor_pressure(300, allow_extrapolation=True)
        self.assertEqual(result.warnings, ("vapor_pressure extrapolated outside 83.65..191.03 K",))
        with self.assertRaises(ValidationError):
            ideal_gas_density(16.04246, 0, 101325)
        with self.assertRaises(ValidationError):
            methane.heat_capacity(0, allow_extrapolation=True)

    def test_correlation_data_retains_provenance_and_matches_vendored_chemsep(self):
        methane = self.correlations["Methane"]
        self.assertEqual(methane.provenance.source_revision, "9.0.4")
        source = ROOT / methane.provenance.source
        xml = ElementTree.parse(source).getroot()
        for compound in xml.iter("compound"):
            name = compound.find("CompoundID")
            if name is None or name.attrib.get("value") not in self.correlations:
                continue
            loaded = self.correlations[name.attrib["value"]]
            for xml_name, correlation in (("IdealGasHeatCapacityCp", loaded.heat_capacity_correlation), ("VaporPressure", loaded.vapor_pressure_correlation)):
                values = {node.tag: node.attrib["value"] for node in compound.find(xml_name)}
                self.assertEqual(correlation.equation, int(values["eqno"]))
                self.assertEqual(correlation.coefficients, tuple(float(values[key]) for key in "ABCDE"))
                self.assertEqual((correlation.minimum_k, correlation.maximum_k), (float(values["Tmin"]), float(values["Tmax"])))

    def test_ideal_properties_match_captured_dwsim_results(self):
        captured = json.loads((ROOT / "tests/golden/compound-catalog.json").read_text(encoding="utf-8-sig"))["inputs"]["compounds"]
        molecular_weights = {record["id"]: record["molecular_weight"]["value"] for record in captured}
        for record in captured:
            reference = record["ideal_reference"]
            correlation = self.correlations[record["id"]]
            temperature = reference["heat_capacity_temperature"]["value"]
            expected_cp = correlation.heat_capacity(temperature).value / molecular_weights[record["id"]] / 1000
            self.assertTrue(math.isclose(reference["heat_capacity"]["value"], expected_cp, rel_tol=1e-10))
            vapor_temperature = reference["vapor_pressure_temperature"]["value"]
            self.assertTrue(math.isclose(reference["vapor_pressure"]["value"], correlation.vapor_pressure(vapor_temperature).value, rel_tol=1e-10))


if __name__ == "__main__":
    unittest.main()
