import copy
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim.compounds import load_compounds, load_pr_interactions
from mesim.errors import ValidationError
from mesim.reactions import load_reaction_data
from mesim.thermo.ideal import load_correlations
from mesim.unitops.reactors import (
    continuous_stirred_tank_reactor,
    conversion_reactor,
    equilibrium_reactor,
    gibbs_reactor,
    plug_flow_reactor,
)


ROOT = Path(__file__).parents[1]
DWSIM_REACTOR_PHASE_AGGREGATE_REL_TOL = 2e-7
DWSIM_EQUILIBRIUM_COMPONENT_REL_TOL = 2e-5
DWSIM_EQUILIBRIUM_EXTENT_REL_TOL = 2e-6
DWSIM_EQUILIBRIUM_DUTY_REL_TOL = 2e-6
DWSIM_GIBBS_COMPONENT_REL_TOL = 1.2e-3
DWSIM_GIBBS_DUTY_REL_TOL = 2e-4
DWSIM_CSTR_RATE_REL_TOL = 3e-5
DWSIM_CSTR_STREAM_REL_TOL = 1.1e-3
DWSIM_PFR_REL_TOL = 2e-4


class ReactionDataTest(unittest.TestCase):
    def test_isomerisation_data_is_explicit_balanced_and_source_backed(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        reaction = data.reactions[0]
        thermo = {record.compound_id: record for record in data.thermochemistry}
        golden = json.loads((ROOT / "tests/golden/u4-conversion-reactor-isomerization.json").read_text(encoding="utf-8-sig"))
        captured = {record["id"]: record for record in golden["inputs"]["compounds"]}

        self.assertEqual(dict(thermo["N-butane"].elements), {"C": 4.0, "H": 10.0})
        self.assertEqual(dict(thermo["Isobutane"].elements), {"C": 4.0, "H": 10.0})
        self.assertEqual(dict(reaction.stoichiometry), {"N-butane": -1.0, "Isobutane": 1.0})
        self.assertEqual(reaction.reaction_heat_j_per_kmol, -9_200_000.0)
        self.assertEqual(reaction.conversion_fraction, 0.33)
        self.assertEqual(data.provenance.source_revision, "9.0.5.0")
        for compound_id in ("N-butane", "Isobutane"):
            record = thermo[compound_id]
            reference = captured[compound_id]
            molecular_weight = reference["molecular_weight"]["value"]
            self.assertEqual(dict(record.elements), reference["elements"])
            self.assertTrue(math.isclose(
                record.ideal_gas_formation_enthalpy_j_per_kmol,
                reference["ideal_gas_formation"]["enthalpy"]["value"] * molecular_weight * 1_000.0,
                rel_tol=1e-12,
            ))
            self.assertTrue(math.isclose(
                record.ideal_gas_formation_gibbs_energy_j_per_kmol,
                reference["ideal_gas_formation"]["gibbs_energy"]["value"] * molecular_weight * 1_000.0,
                rel_tol=1e-12,
            ))

    def test_reaction_loader_rejects_unbalanced_stoichiometry_and_inconsistent_heat(self):
        source = json.loads((ROOT / "data/reactions/v1.json").read_text())
        for mutation in ("stoichiometry", "heat", "formation entropy"):
            broken = json.loads(json.dumps(source))
            if mutation == "stoichiometry":
                broken["reactions"][0]["stoichiometry"]["Isobutane"] = 2.0
            elif mutation == "heat":
                broken["reactions"][0]["reaction_heat"]["value"] = 0.0
            else:
                broken["thermochemistry"][0]["ideal_gas_formation_entropy"]["value"] = 0.0
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "reactions.json"
                path.write_text(json.dumps(broken))
                with self.assertRaises(ValidationError):
                    load_reaction_data(path)

    def test_ethylene_glycol_kinetics_preserve_dwsim_expression_units(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        reaction = next(record for record in data.reactions if record.reaction_type == "kinetic")
        kinetics = reaction.kinetics

        self.assertIsNotNone(kinetics)
        assert kinetics is not None
        self.assertEqual(reaction.phase, "mixture")
        self.assertEqual(dict(reaction.stoichiometry), {
            "Ethylene oxide": -1.0, "Water": -1.0, "Ethylene glycol": 1.0,
        })
        self.assertEqual(reaction.reaction_heat_j_per_kmol, -97_756_000.0)
        self.assertEqual(kinetics.basis, "molar_concentration")
        self.assertEqual(kinetics.concentration_unit, "kmol/m3")
        self.assertEqual(kinetics.rate_unit, "kmol/[m3.h]")
        self.assertEqual(kinetics.forward.pre_exponential_factor, 0.005)
        self.assertEqual(kinetics.forward.activation_energy_j_per_mol, 0.0)
        self.assertEqual(dict(kinetics.forward.orders)["Ethylene oxide"], 1.0)

        golden = json.loads(
            (ROOT / "tests/golden/u5-cstr-ethylene-glycol-raoult.json").read_text(encoding="utf-8-sig")
        )
        captured = {record["id"]: record for record in golden["inputs"]["compounds"]}
        thermochemistry = {record.compound_id: record for record in data.thermochemistry}
        for compound_id in ("Ethylene oxide", "Ethylene glycol"):
            record = thermochemistry[compound_id]
            reference = captured[compound_id]
            molecular_weight = reference["molecular_weight"]["value"]
            self.assertEqual(dict(record.elements), reference["elements"])
            self.assertTrue(math.isclose(
                record.ideal_gas_formation_enthalpy_j_per_kmol,
                reference["ideal_gas_formation"]["enthalpy"]["value"] * molecular_weight * 1_000.0,
                rel_tol=1e-12,
            ))

        with ZipFile(ROOT / "tests/u5-cstr-ethylene-glycol-raoult.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        saved = root.find(".//Reaction[ReactionType='Kinetic']")
        self.assertIsNotNone(saved)
        assert saved is not None
        self.assertEqual(saved.findtext("ConcUnit"), kinetics.concentration_unit)
        self.assertEqual(saved.findtext("VelUnit"), kinetics.rate_unit)
        self.assertEqual(saved.findtext("E_Forward_Unit"), kinetics.forward.activation_energy_unit)
        self.assertEqual(float(saved.findtext("A_Forward")), kinetics.forward.pre_exponential_factor)
        self.assertEqual(float(saved.findtext("E_Forward")), kinetics.forward.activation_energy_j_per_mol)


class ConversionReactorTest(unittest.TestCase):
    def test_conversion_reactor_matches_captured_dwsim_isomerisation(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        reaction = data.reactions[0]
        golden = json.loads((ROOT / "tests/golden/u4-conversion-reactor-isomerization.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in objects[tag]["properties"]}
            for tag in ("R4", "Reactor Output1", "Reactor Output 2", "Conversion Reactor")
        }
        compounds = ("Propane", "Isobutane", "N-butane", "Isopentane")
        inlet = tuple(
            (compound, properties["R4"][f"PROP_MS_104/{compound}"] / 1_000.0)
            for compound in compounds
        )
        expected_outlet = {
            compound: math.fsum(
                properties[tag][f"PROP_MS_104/{compound}"] / 1_000.0
                for tag in ("Reactor Output1", "Reactor Output 2")
            )
            for compound in compounds
        }
        result = conversion_reactor(inlet, reaction)
        actual_outlet = dict(result.outlet_component_flows_kmol_s)

        self.assertEqual(properties["Conversion Reactor"]["isomerisation: Extent"], 33.0)
        self.assertEqual(properties["Conversion Reactor"]["N-butane: Conversion"], 33.0)
        self.assertEqual(result.extent_kmol_s, dict(inlet)["N-butane"] * 0.33)
        self.assertEqual(result.reaction_heat_w, result.extent_kmol_s * -9_200_000.0)
        for compound in compounds:
            self.assertTrue(math.isclose(
                actual_outlet[compound], expected_outlet[compound],
                rel_tol=DWSIM_REACTOR_PHASE_AGGREGATE_REL_TOL,
            ))
        inlet_total = math.fsum(flow for _, flow in inlet)
        expected_total = math.fsum(expected_outlet.values())
        self.assertTrue(math.isclose(result.total_molar_flow_kmol_s, inlet_total, rel_tol=0.0, abs_tol=1e-15))
        self.assertTrue(math.isclose(expected_total, inlet_total, rel_tol=1e-9))

        with ZipFile(ROOT / "tests/u4-conversion-reactor-isomerization.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        saved = root.find("./Reactions/Reaction")
        self.assertEqual(saved.findtext("ReactionType"), "Conversion")
        self.assertEqual(saved.findtext("BaseReactant"), reaction.base_reactant)
        self.assertEqual(float(saved.findtext("Expression")) / 100.0, reaction.conversion_fraction)
        self.assertTrue(math.isclose(float(saved.findtext("ReactionHeat")) * 1_000.0, reaction.reaction_heat_j_per_kmol, rel_tol=1e-12))

    def test_conversion_reactor_rejects_invalid_conversion_and_missing_base_flow(self):
        reaction = load_reaction_data(ROOT / "data/reactions/v1.json").reactions[0]
        with self.assertRaises(ValidationError):
            conversion_reactor((("N-butane", 1.0),), reaction, 1.01)
        with self.assertRaises(ValidationError):
            conversion_reactor((("Isobutane", 1.0),), reaction)


class EquilibriumReactorTest(unittest.TestCase):
    def test_steam_reforming_matches_captured_dwsim_equilibrium_reactor(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        reactions = tuple(reaction for reaction in data.reactions if reaction.reaction_type == "equilibrium")
        catalog = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}
        correlation_catalog = {
            record.compound_id: record for record in load_correlations(ROOT / "data/correlations/ideal-v1.json")
        }
        interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        golden = json.loads(
            (ROOT / "tests/golden/u4-equilibrium-reactor-steam-reforming-pr-eos.json").read_text(encoding="utf-8-sig")
        )
        repeat = json.loads(
            (ROOT / "tests/golden/u4-equilibrium-reactor-steam-reforming-pr-eos-repeat.json").read_text(encoding="utf-8-sig")
        )
        normalized_golden, normalized_repeat = copy.deepcopy(golden), copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in objects[tag]["properties"]}
            for tag in ("4", "5", "e2", "RE-001")
        }
        compound_ids = ("Methane", "Water", "Carbon monoxide", "Carbon dioxide", "Hydrogen")
        inlet = tuple(
            (compound_id, properties["4"][f"PROP_MS_104/{compound_id}"] / 1_000.0)
            for compound_id in compound_ids
        )
        result = equilibrium_reactor(
            inlet,
            reactions,
            data.thermochemistry,
            tuple(catalog[compound_id] for compound_id in compound_ids),
            tuple(correlation_catalog[compound_id] for compound_id in compound_ids),
            interactions,
            properties["4"]["PROP_MS_0"],
            properties["4"]["PROP_MS_1"],
        )

        expected_outlet = {
            compound_id: properties["5"][f"PROP_MS_104/{compound_id}"] / 1_000.0
            for compound_id in compound_ids
        }
        for compound_id, actual in result.outlet_component_flows_kmol_s:
            self.assertTrue(math.isclose(
                actual, expected_outlet[compound_id],
                rel_tol=DWSIM_EQUILIBRIUM_COMPONENT_REL_TOL, abs_tol=1.0e-10,
            ))
        for reaction, (_, extent) in zip(reactions, result.extents_kmol_s):
            self.assertTrue(math.isclose(
                extent, properties["RE-001"][f"{reaction.name}: Extent"] / 1_000.0,
                rel_tol=DWSIM_EQUILIBRIUM_EXTENT_REL_TOL,
            ))
        conversions = dict(result.component_conversions)
        self.assertTrue(math.isclose(
            conversions["Methane"] * 100.0, properties["RE-001"]["Methane: Conversion"],
            rel_tol=DWSIM_EQUILIBRIUM_COMPONENT_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            conversions["Water"] * 100.0, properties["RE-001"]["Water: Conversion"],
            rel_tol=DWSIM_EQUILIBRIUM_COMPONENT_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            result.isothermal_duty_w / 1_000.0, properties["e2"]["PROP_ES_0"],
            rel_tol=DWSIM_EQUILIBRIUM_DUTY_REL_TOL,
        ))
        self.assertLess(max(abs(value) for _, value in result.equilibrium_log_residuals), 1.0e-9)

    def test_equilibrium_reactor_rejects_a_conversion_reaction(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        with self.assertRaises(ValidationError):
            equilibrium_reactor(
                (("N-butane", 1.0), ("Isobutane", 0.1)),
                (data.reactions[0],), data.thermochemistry, (), (),
                load_pr_interactions(ROOT / "data/interactions/pr-v1.json"), 400.0, 101325.0,
            )


class CSTRReactorTest(unittest.TestCase):
    def test_ethylene_glycol_cstr_matches_repeatable_dwsim_case(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        reaction = next(record for record in data.reactions if record.reaction_type == "kinetic")
        golden = json.loads(
            (ROOT / "tests/golden/u5-cstr-ethylene-glycol-raoult.json").read_text(encoding="utf-8-sig")
        )
        repeat = json.loads(
            (ROOT / "tests/golden/u5-cstr-ethylene-glycol-raoult-repeat.json").read_text(encoding="utf-8-sig")
        )
        normalized_golden, normalized_repeat = copy.deepcopy(golden), copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        self.assertFalse(any(
            property_record["read_error"]
            for object_record in golden["outputs"]["objects_after"]
            for property_record in object_record["properties"]
        ))

        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in objects[tag]["properties"]}
            for tag in ("5", "6", "E2", "CSTR-1")
        }
        compound_ids = ("Ethylene oxide", "Water", "Ethylene glycol")
        inlet = tuple(
            (compound_id, properties["5"][f"PROP_MS_104/{compound_id}"] / 1_000.0)
            for compound_id in compound_ids
        )
        result = continuous_stirred_tank_reactor(
            inlet,
            reaction,
            properties["5"]["PROP_MS_0"],
            properties["CSTR-1"]["PROP_CS_2"],
            properties["6"]["PROP_MS_4"],
        )

        expected_outlet = {
            compound_id: properties["6"][f"PROP_MS_104/{compound_id}"] / 1_000.0
            for compound_id in compound_ids
        }
        for compound_id, actual in result.outlet_component_flows_kmol_s:
            self.assertTrue(math.isclose(
                actual, expected_outlet[compound_id],
                rel_tol=DWSIM_CSTR_STREAM_REL_TOL, abs_tol=1.0e-12,
            ))
        self.assertTrue(math.isclose(
            result.extent_kmol_s * 1_000.0,
            properties["CSTR-1"]["Ethylene Glycol Production: Extent"],
            rel_tol=DWSIM_CSTR_RATE_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            result.reaction_rate_kmol_m3_s * 1_000.0,
            properties["CSTR-1"]["Ethylene Glycol Production: Rate"],
            rel_tol=DWSIM_CSTR_RATE_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            result.reference_reaction_heat_w / 1_000.0,
            properties["CSTR-1"]["Ethylene Glycol Production: Heat"],
            rel_tol=DWSIM_CSTR_RATE_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            dict(result.component_conversions)["Ethylene oxide"] * 100.0,
            properties["CSTR-1"]["Ethylene oxide: Conversion"],
            rel_tol=DWSIM_CSTR_STREAM_REL_TOL,
        ))
        self.assertLess(abs(result.material_rate_residual_kmol_s), 1.0e-15)
        self.assertTrue(math.isclose(
            result.total_molar_flow_kmol_s,
            math.fsum(flow for _, flow in inlet) - result.extent_kmol_s,
            rel_tol=0.0, abs_tol=1.0e-15,
        ))

    def test_cstr_rejects_nonkinetic_reaction_and_invalid_supplied_volume(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        kinetic = next(record for record in data.reactions if record.reaction_type == "kinetic")
        with self.assertRaises(ValidationError):
            continuous_stirred_tank_reactor(
                (("N-butane", 1.0), ("Isobutane", 0.0)), data.reactions[0], 350.0, 1.0, 0.001,
            )
        with self.assertRaises(ValidationError):
            continuous_stirred_tank_reactor(
                (("Ethylene oxide", 1.0), ("Water", 1.0), ("Ethylene glycol", 0.0)),
                kinetic, 350.0, 1.0, 0.0,
            )


class PFRReactorTest(unittest.TestCase):
    def test_ethylene_glycol_pfr_matches_repeatable_dwsim_case(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        reaction = next(record for record in data.reactions if record.reaction_type == "kinetic")
        golden = json.loads(
            (ROOT / "tests/golden/u5-pfr-ethylene-glycol-raoult.json").read_text(encoding="utf-8-sig")
        )
        repeat = json.loads(
            (ROOT / "tests/golden/u5-pfr-ethylene-glycol-raoult-repeat.json").read_text(encoding="utf-8-sig")
        )
        normalized_golden, normalized_repeat = copy.deepcopy(golden), copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        self.assertFalse(any(
            property_record["read_error"]
            for object_record in golden["outputs"]["objects_after"]
            for property_record in object_record["properties"]
        ))

        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in objects[tag]["properties"]}
            for tag in ("Feed", "Product", "Energy", "PFR-001")
        }
        compound_ids = ("Ethylene oxide", "Water", "Ethylene glycol")
        inlet = tuple(
            (compound_id, properties["Feed"][f"PROP_MS_104/{compound_id}"] / 1_000.0)
            for compound_id in compound_ids
        )
        result = plug_flow_reactor(
            inlet,
            reaction,
            properties["Feed"]["PROP_MS_0"],
            properties["Product"]["PROP_MS_0"],
            properties["PFR-001"]["PROP_PF_2"],
            properties["Feed"]["PROP_MS_4"],
            properties["Product"]["PROP_MS_4"],
        )

        expected_outlet = {
            compound_id: properties["Product"][f"PROP_MS_104/{compound_id}"] / 1_000.0
            for compound_id in compound_ids
        }
        for compound_id, actual in result.outlet_component_flows_kmol_s:
            self.assertTrue(math.isclose(
                actual, expected_outlet[compound_id],
                rel_tol=DWSIM_PFR_REL_TOL, abs_tol=1.0e-12,
            ))
        # DWSIM labels these two properties as mol units, but their numeric values
        # close the material and heat balances only when read as kmol/s and kmol/m3/s.
        self.assertTrue(math.isclose(
            result.extent_kmol_s,
            properties["PFR-001"]["Ethylene Glycol Production: Reaction Extent"],
            rel_tol=DWSIM_PFR_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            result.average_reaction_rate_kmol_m3_s,
            properties["PFR-001"]["Ethylene Glycol Production: Reaction Rate"],
            rel_tol=DWSIM_PFR_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            result.reference_reaction_heat_w / 1_000.0,
            properties["PFR-001"]["Ethylene Glycol Production: Reaction Heat"],
            rel_tol=DWSIM_PFR_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            dict(result.component_conversions)["Ethylene oxide"] * 100.0,
            properties["PFR-001"]["Ethylene oxide: Conversion"],
            rel_tol=DWSIM_PFR_REL_TOL,
        ))
        self.assertLess(result.material_balance_residual_kmol_s, 1.0e-15)
        self.assertTrue(math.isclose(
            result.total_molar_flow_kmol_s,
            math.fsum(flow for _, flow in inlet) - result.extent_kmol_s,
            rel_tol=0.0, abs_tol=1.0e-15,
        ))

    def test_pfr_rejects_nonkinetic_reaction_and_invalid_profile(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        kinetic = next(record for record in data.reactions if record.reaction_type == "kinetic")
        inlet = (("Ethylene oxide", 1.0), ("Water", 1.0), ("Ethylene glycol", 0.0))
        with self.assertRaises(ValidationError):
            plug_flow_reactor(inlet, data.reactions[0], 350.0, 350.0, 1.0, 0.001, 0.001)
        with self.assertRaises(ValidationError):
            plug_flow_reactor(inlet, kinetic, 350.0, 350.0, 1.0, 0.001, 0.0)


class GibbsReactorTest(unittest.TestCase):
    def test_steam_reforming_minimization_matches_dwsim_and_closes_elements(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        catalog = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}
        correlation_catalog = {
            record.compound_id: record for record in load_correlations(ROOT / "data/correlations/ideal-v1.json")
        }
        golden = json.loads(
            (ROOT / "tests/golden/u4-gibbs-reactor-steam-reforming-pr-eos.json").read_text(encoding="utf-8-sig")
        )
        repeat = json.loads(
            (ROOT / "tests/golden/u4-gibbs-reactor-steam-reforming-pr-eos-repeat.json").read_text(encoding="utf-8-sig")
        )
        normalized_golden, normalized_repeat = copy.deepcopy(golden), copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in objects[tag]["properties"]}
            for tag in ("1", "2", "e1", "RG-000")
        }
        compound_ids = ("Methane", "Water", "Carbon monoxide", "Carbon dioxide", "Hydrogen")
        inlet = tuple(
            (compound_id, properties["1"][f"PROP_MS_104/{compound_id}"] / 1_000.0)
            for compound_id in compound_ids
        )
        result = gibbs_reactor(
            inlet,
            data.thermochemistry,
            tuple(catalog[compound_id] for compound_id in compound_ids),
            tuple(correlation_catalog[compound_id] for compound_id in compound_ids),
            load_pr_interactions(ROOT / "data/interactions/pr-v1.json"),
            properties["1"]["PROP_MS_0"],
            properties["1"]["PROP_MS_1"],
        )

        expected_outlet = {
            compound_id: properties["2"][f"PROP_MS_104/{compound_id}"] / 1_000.0
            for compound_id in compound_ids
        }
        for compound_id, actual in result.outlet_component_flows_kmol_s:
            self.assertTrue(math.isclose(
                actual, expected_outlet[compound_id],
                rel_tol=DWSIM_GIBBS_COMPONENT_REL_TOL, abs_tol=1.0e-10,
            ))
        conversions = dict(result.component_conversions)
        self.assertTrue(math.isclose(
            conversions["Methane"] * 100.0, properties["RG-000"]["Methane: Conversion"],
            rel_tol=DWSIM_GIBBS_COMPONENT_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            conversions["Water"] * 100.0, properties["RG-000"]["Water: Conversion"],
            rel_tol=DWSIM_GIBBS_COMPONENT_REL_TOL,
        ))
        self.assertTrue(math.isclose(
            result.isothermal_duty_w / 1_000.0, properties["e1"]["PROP_ES_0"],
            rel_tol=DWSIM_GIBBS_DUTY_REL_TOL,
        ))
        self.assertLess(max(abs(value) for _, value in result.element_balance_residuals_kmol_s), 1.0e-11)
        self.assertLess(result.stationarity_residual, 1.0e-9)
        self.assertLess(result.final_gibbs_energy_w, result.initial_gibbs_energy_w)

        captured_oxygen = (
            expected_outlet["Water"] + expected_outlet["Carbon monoxide"]
            + 2.0 * expected_outlet["Carbon dioxide"]
        )
        inlet_oxygen = dict(inlet)["Water"]
        self.assertGreater(abs(captured_oxygen - inlet_oxygen), 1.0e-7)

    def test_gibbs_reactor_rejects_zero_flow_and_missing_data(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        with self.assertRaises(ValidationError):
            gibbs_reactor((("Methane", 0.0), ("Water", 0.0)), data.thermochemistry, (), (), interactions, 1000.0, 101325.0)


if __name__ == "__main__":
    unittest.main()
