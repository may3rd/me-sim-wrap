import unittest
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.compounds import load_compounds, load_pr_interactions
from mesim.errors import ValidationError
from mesim.thermo.activity import (
    load_nrtl_vle_data,
    nrtl_bubble_pressure,
    nrtl_equilibrium_ratios,
    nrtl_phase_enthalpies,
)
from mesim.thermo.flash import flash_enthalpy, tp_flash
from mesim.thermo.ideal import load_correlations
from mesim.thermo.systems import (
    IDEAL_RAOULT,
    NRTL_ACETONE_METHANOL,
    PENG_ROBINSON_CLASSIC,
    PENG_ROBINSON_1978,
    PENG_ROBINSON_LEE_KESLER,
    PENG_ROBINSON_STRYJEK_VERA_2_MARGULES,
    SOAVE_REDLICH_KWONG,
    THERMO_SYSTEM_CONSTRUCTORS,
    IdealRaoultSystem,
    NRTLSystem,
    PengRobinsonSystem,
    ThermodynamicSystem,
    create_thermo_system,
)


class ThermodynamicSystemTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        catalog = {
            compound.id: compound
            for compound in load_compounds(ROOT / "data/compounds/v1.json")
        }
        cls.compounds = (catalog["Methane"], catalog["Ethane"])
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        cls.correlations = load_correlations(ROOT / "data/correlations/ideal-v1.json")
        cls.nrtl_data = load_nrtl_vle_data(
            ROOT / "data/correlations/nrtl-acetone-methanol-v1.json"
        )

    def test_registry_has_stable_non_plugin_model_ids(self):
        self.assertEqual(
            set(THERMO_SYSTEM_CONSTRUCTORS),
            {
                PENG_ROBINSON_CLASSIC,
                NRTL_ACETONE_METHANOL,
                IDEAL_RAOULT,
                SOAVE_REDLICH_KWONG,
                PENG_ROBINSON_1978,
                PENG_ROBINSON_LEE_KESLER,
                PENG_ROBINSON_STRYJEK_VERA_2_MARGULES,
            },
        )
        with self.assertRaises(TypeError):
            THERMO_SYSTEM_CONSTRUCTORS["runtime-plugin"] = object
        with self.assertRaises(ValidationError):
            create_thermo_system("unknown-model")
        with self.assertRaises(ValidationError):
            create_thermo_system(PENG_ROBINSON_CLASSIC, compounds=self.compounds)

    def test_pr_system_preserves_flash_and_caloric_results(self):
        system = create_thermo_system(
            PENG_ROBINSON_CLASSIC,
            compounds=self.compounds,
            interactions=self.interactions,
            correlations=self.correlations,
        )
        self.assertIsInstance(system, PengRobinsonSystem)
        self.assertIsInstance(system, ThermodynamicSystem)
        self.assertEqual(system.compound_ids, ("Methane", "Ethane"))

        composition = (0.7, 0.3)
        direct = tp_flash(
            self.compounds, composition, self.interactions, 180.0, 500_000.0
        )
        extracted = system.tp_flash(composition, 180.0, 500_000.0)
        self.assertEqual(extracted, direct)
        expected_enthalpy = flash_enthalpy(
            self.compounds, self.correlations, direct
        )
        self.assertEqual(system.enthalpy(extracted), expected_enthalpy)
        round_trip = system.ph_flash(
            composition,
            500_000.0,
            expected_enthalpy,
            (140.0, 210.0),
        )
        self.assertTrue(round_trip.report.converged)
        self.assertAlmostEqual(round_trip.temperature_k, 180.0, places=6)

    def test_nrtl_system_preserves_stage_equilibrium_and_calorics(self):
        ids = ("Methanol", "Acetone")
        system = create_thermo_system(
            NRTL_ACETONE_METHANOL,
            data=self.nrtl_data,
            compound_ids=ids,
        )
        self.assertIsInstance(system, NRTLSystem)
        self.assertIsInstance(system, ThermodynamicSystem)
        composition = (0.30481315, 0.69518685)
        temperature_k = 388.288289
        pressure_pa = 607_950.0

        expected_ratios = nrtl_equilibrium_ratios(
            self.nrtl_data, ids, composition, temperature_k, pressure_pa
        )
        self.assertEqual(
            system.equilibrium_ratios(composition, temperature_k, pressure_pa),
            expected_ratios,
        )
        expected_bubble = nrtl_bubble_pressure(
            self.nrtl_data, ids, composition, temperature_k
        )
        bubble = system.bubble_pressure(composition, temperature_k)
        self.assertEqual(bubble, expected_bubble)
        expected_enthalpies = nrtl_phase_enthalpies(
            self.nrtl_data,
            ids,
            composition,
            bubble.vapor_composition,
            temperature_k,
            pressure_pa,
        )
        self.assertEqual(
            system.phase_enthalpies(
                composition,
                bubble.vapor_composition,
                temperature_k,
                pressure_pa,
            ),
            expected_enthalpies,
        )

    def test_ideal_raoult_system_preserves_direct_equilibrium(self):
        catalog = {record.compound_id: record for record in self.correlations}
        correlations = (catalog["Methane"], catalog["Ethane"])
        system = create_thermo_system(IDEAL_RAOULT, correlations=correlations)
        self.assertIsInstance(system, IdealRaoultSystem)
        self.assertIsInstance(system, ThermodynamicSystem)
        self.assertEqual(system.compound_ids, ("Methane", "Ethane"))

        temperature_k = 180.0
        pressure_pa = 500_000.0
        self.assertEqual(
            system.fugacity_coefficients(temperature_k, pressure_pa, "vapor"),
            (1.0, 1.0),
        )
        bubble = system.bubble_pressure((0.7, 0.3), temperature_k)
        self.assertTrue(bubble.report.converged)
        self.assertEqual(
            bubble.equilibrium_ratios,
            system.equilibrium_ratios(temperature_k, bubble.pressure_pa),
        )

    def test_systems_reject_incomplete_or_mismatched_domains(self):
        without_ethane = tuple(
            record for record in self.correlations if record.compound_id != "Ethane"
        )
        with self.assertRaises(ValidationError):
            PengRobinsonSystem(
                self.compounds, self.interactions, without_ethane
            )
        with self.assertRaises(ValidationError):
            NRTLSystem(self.nrtl_data, ("Methanol", "Unknown"))


if __name__ == "__main__":
    unittest.main()
