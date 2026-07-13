import json
import math
import sys
import unittest
from xml.etree import ElementTree
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import OutOfRangeError, ValidationError
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.thermo.ideal import ideal_gas_density, load_correlations
from mesim.thermo.peng_robinson import PengRobinson, PengRobinsonMixture


ROOT = Path(__file__).parents[1]
DWSIM_PR_FUGACITY_REL_TOL = 2e-6
DWSIM_PR_CALORIC_DELTA_REL_TOL = 2e-4


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
        extreme_density = ideal_gas_density(1e200, 1e100, 1e200).value
        self.assertTrue(math.isfinite(extreme_density))
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


class PengRobinsonPureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compounds = load_compounds(ROOT / "data/compounds/v1.json")
        cls.methane = PengRobinson(next(c for c in compounds if c.id == "Methane"))

    def test_methane_parameters_match_independent_pr_equations(self):
        parameters = self.methane.parameters(300)
        self.assertTrue(math.isclose(parameters.kappa, 0.39157219968, rel_tol=1e-13))
        self.assertTrue(math.isclose(parameters.alpha, 0.8104699865610172, rel_tol=1e-13))
        self.assertTrue(math.isclose(parameters.a_pa_m6_per_kmol2, 202278.44274980662, rel_tol=1e-13))
        self.assertTrue(math.isclose(parameters.b_m3_per_kmol, 0.026802920402019762, rel_tol=1e-13))

    def test_cubic_roots_are_real_physical_and_phase_selected(self):
        roots = self.methane.roots(150, 1_000_000)
        expected = (0.03312399834715224, 0.12030482166960407, 0.8250801775914196)
        self.assertEqual(len(roots), 3)
        self.assertTrue(all(math.isclose(actual, wanted, rel_tol=1e-12) for actual, wanted in zip(roots, expected)))
        self.assertEqual(self.methane.state(150, 1_000_000, "liquid").compressibility, roots[0])
        self.assertEqual(self.methane.state(150, 1_000_000, "vapor").compressibility, roots[-1])
        self.assertEqual(self.methane.stable_state(150, 800_000).phase, "vapor")
        self.assertEqual(self.methane.stable_state(150, 1_200_000).phase, "liquid")
        self.assertEqual(self.methane.stable_state(300, 1_000_000).phase, "single")

    def test_vapor_state_matches_independent_pr_properties(self):
        state = self.methane.state(300, 1_000_000, "vapor")
        self.assertTrue(math.isclose(state.compressibility, 0.9785896699519482, rel_tol=1e-12))
        self.assertTrue(math.isclose(state.fugacity_coefficient, 0.9786405888609191, rel_tol=1e-12))
        self.assertTrue(math.isclose(state.density_kg_per_m3, 6.5722624579836175, rel_tol=1e-12))
        self.assertTrue(math.isclose(state.departure_enthalpy_j_per_kmol, -180117.9009712916, rel_tol=1e-12))
        self.assertTrue(math.isclose(state.departure_entropy_j_per_kmol_k, -420.87689978670664, rel_tol=1e-12))

    def test_invalid_pure_states_are_rejected(self):
        for temperature, pressure in (
            (0, 1), (300, 0), (math.nan, 1), (300, math.inf),
            (300, 5e-324), (300, 1e12), (1e-300, 1e300), (1e300, 1e-300),
            (1e100, 1e200), (1e-200, 1e-200),
        ):
            with self.assertRaises(ValidationError):
                self.methane.state(temperature, pressure, "vapor")
        with self.assertRaises(ValidationError):
            self.methane.state(300, 1_000_000, "unknown")

    def test_pure_fugacity_matches_dwsim_pr_vapor_and_liquid_cases(self):
        golden = json.loads((ROOT / "tests/golden/u0-pr-c1-c5.json").read_text(encoding="utf-8-sig"))
        streams = {record["tag"]: record for record in golden["outputs"]["objects_after"]}

        def property_value(tag, name):
            properties = {record["property"]: record["value"]["value"] for record in streams[tag]["properties"]}
            return properties[name]

        methane_reference = property_value("CH4-feed", "Fugacity Coefficient, Vapor Phase / Methane")
        self.assertTrue(math.isclose(self.methane.state(300, 1_000_000, "vapor").fugacity_coefficient, methane_reference, rel_tol=1e-10))

        compounds = load_compounds(ROOT / "data/compounds/v1.json")
        pentane = PengRobinson(next(c for c in compounds if c.id == "N-pentane"))
        pentane_reference = property_value("C5-feed", "Fugacity Coefficient, Liquid Phase 1 / N-pentane")
        self.assertTrue(math.isclose(pentane.state(300, 1_000_000, "liquid").fugacity_coefficient, pentane_reference, rel_tol=2e-6))


class PengRobinsonMixtureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compounds = load_compounds(ROOT / "data/compounds/v1.json")
        selected = tuple(c for c in compounds if c.id in {"Methane", "Ethane"})
        interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        cls.mixture = PengRobinsonMixture(selected, (0.7, 0.3), interactions)

    def test_mixing_rules_match_independent_equations(self):
        parameters = self.mixture.parameters(250)
        self.assertTrue(math.isclose(parameters.a_pa_m6_per_kmol2, 330869.45686457393, rel_tol=1e-12))
        self.assertTrue(math.isclose(parameters.b_m3_per_kmol, 0.030923428538033274, rel_tol=1e-12))
        self.assertTrue(math.isclose(parameters.da_dtemperature, -616.2303397589949, rel_tol=1e-12))

    def test_mixture_state_matches_independent_pr_equations(self):
        state = self.mixture.state(250, 5_000_000, "vapor")
        self.assertTrue(math.isclose(state.compressibility, 0.6411750398839282, rel_tol=1e-12))
        self.assertTrue(math.isclose(state.fugacity_coefficients[0], 0.8525219448372683, rel_tol=1e-12))
        self.assertTrue(math.isclose(state.fugacity_coefficients[1], 0.4829182900524611, rel_tol=1e-12))
        self.assertTrue(math.isclose(state.density_kg_per_m3, 75.9719962978659, rel_tol=1e-12))
        self.assertTrue(math.isclose(state.departure_enthalpy_j_per_kmol, -2387903.066060099, rel_tol=1e-12))
        self.assertTrue(math.isclose(state.departure_entropy_j_per_kmol_k, -6807.325988454073, rel_tol=1e-12))
        self.assertEqual(self.mixture.stable_state(250, 5_000_000).phase, "single")

    def test_composition_is_validated_at_construction(self):
        compounds = tuple(model.compound for model in self.mixture.components)
        interactions = self.mixture.interactions
        for fractions in ((0.7,), (0.7, 0.4), (-0.1, 1.1), (math.nan, math.nan)):
            with self.assertRaises(ValidationError):
                PengRobinsonMixture(compounds, fractions, interactions)
        for temperature, pressure in (
            (300, 5e-324), (300, 1e12), (1e-300, 1e300), (1e300, 1e-300),
            (1e100, 1e200), (1e-200, 1e-200),
        ):
            with self.assertRaises(ValidationError):
                self.mixture.state(temperature, pressure, "vapor")

    def test_pr_properties_match_captured_dwsim_reference_states(self):
        golden = json.loads((ROOT / "tests/golden/pr-t1.json").read_text(encoding="utf-8-sig"))
        self.assertEqual(golden["outputs"]["solve"], {"errors": [], "success": True, "executed": True})
        streams = {
            record["tag"]: {item["property"]: item["value"]["value"] for item in record["properties"]}
            for record in golden["outputs"]["objects_after"]
        }
        pure_cases = {
            "PR-V-METHANE": ("Methane", "vapor"),
            "PR-V-ETHANE": ("Ethane", "vapor"),
            "PR-V-PROPANE": ("Propane", "vapor"),
            "PR-V-NBUTANE": ("N-butane", "vapor"),
            "PR-V-NPENTANE": ("N-pentane", "vapor"),
            "PR-L-METHANE": ("Methane", "liquid"),
            "PR-L-ETHANE": ("Ethane", "liquid"),
            "PR-L-PROPANE": ("Propane", "liquid"),
            "PR-L-NBUTANE": ("N-butane", "liquid"),
            "PR-L-NPENTANE": ("N-pentane", "liquid"),
            "PR-3ROOT-METHANE": ("Methane", "vapor"),
            "PR-NC-METHANE": ("Methane", "liquid"),
        }
        compounds = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}

        for tag, (compound_id, phase) in pure_cases.items():
            with self.subTest(tag=tag):
                properties = streams[tag]
                model = PengRobinson(compounds[compound_id])
                state = model.state(properties["PROP_MS_0"], properties["PROP_MS_1"], phase)
                phase_name = "Vapor Phase" if phase == "vapor" else "Liquid Phase 1"
                reference = properties[f"Fugacity Coefficient, {phase_name} / {compound_id}"]
                self.assertTrue(math.isclose(state.fugacity_coefficient, reference, rel_tol=DWSIM_PR_FUGACITY_REL_TOL))

        for tag, expected_phase in (("PR-3ROOT-METHANE", "vapor"), ("PR-NC-METHANE", "liquid")):
            properties = streams[tag]
            model = PengRobinson(compounds["Methane"])
            self.assertEqual(model.stable_state(properties["PROP_MS_0"], properties["PROP_MS_1"]).phase, expected_phase)

        interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        mixture = PengRobinsonMixture((compounds["Methane"], compounds["Ethane"]), (0.7, 0.3), interactions)
        mixture_properties = streams["PR-MIX-ME-C2"]
        mixture_state = mixture.state(mixture_properties["PROP_MS_0"], mixture_properties["PROP_MS_1"], "vapor")
        for compound_id, actual in zip(("Methane", "Ethane"), mixture_state.fugacity_coefficients):
            reference = mixture_properties[f"Fugacity Coefficient, Vapor Phase / {compound_id}"]
            self.assertTrue(math.isclose(actual, reference, rel_tol=DWSIM_PR_FUGACITY_REL_TOL))

        correlations = {record.compound_id: record for record in load_correlations(ROOT / "data/correlations/ideal-v1.json")}
        for suffix, compound_id in (
            ("METHANE", "Methane"),
            ("ETHANE", "Ethane"),
            ("PROPANE", "Propane"),
            ("NBUTANE", "N-butane"),
            ("NPENTANE", "N-pentane"),
        ):
            vapor = streams[f"PR-V-{suffix}"]
            liquid = streams[f"PR-L-{suffix}"]
            model = PengRobinson(compounds[compound_id])
            vapor_state = model.state(vapor["PROP_MS_0"], vapor["PROP_MS_1"], "vapor")
            liquid_state = model.state(liquid["PROP_MS_0"], liquid["PROP_MS_1"], "liquid")
            enthalpy_delta = correlations[compound_id].enthalpy_change(liquid["PROP_MS_0"], vapor["PROP_MS_0"]).value + liquid_state.departure_enthalpy_j_per_kmol - vapor_state.departure_enthalpy_j_per_kmol
            entropy_delta = correlations[compound_id].entropy_change(liquid["PROP_MS_0"], liquid["PROP_MS_1"], vapor["PROP_MS_0"], vapor["PROP_MS_1"]).value + liquid_state.departure_entropy_j_per_kmol_k - vapor_state.departure_entropy_j_per_kmol_k
            self.assertTrue(math.isclose(enthalpy_delta, (liquid["PROP_MS_9"] - vapor["PROP_MS_9"]) * 1_000, rel_tol=DWSIM_PR_CALORIC_DELTA_REL_TOL))
            self.assertTrue(math.isclose(entropy_delta, (liquid["PROP_MS_10"] - vapor["PROP_MS_10"]) * 1_000, rel_tol=DWSIM_PR_CALORIC_DELTA_REL_TOL))


if __name__ == "__main__":
    unittest.main()
