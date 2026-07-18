import json
import math
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.compounds import load_compounds
from mesim.errors import ValidationError
from mesim.thermo.prsv2 import PRSV2Mixture, load_prsv2_data
from mesim.thermo.systems import (
    PENG_ROBINSON_STRYJEK_VERA_2_MARGULES,
    PENG_ROBINSON_STRYJEK_VERA_2_VAN_LAAR,
    PRSV2MargulesSystem,
    PRSV2VanLaarSystem,
    create_thermo_system,
)


class PRSV2Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = {
            compound.id: compound
            for compound in load_compounds(ROOT / "data/compounds/v1.json")
        }
        cls.data_path = ROOT / "data/interactions/prsv2-v1.json"
        cls.data = load_prsv2_data(cls.data_path)

    def test_full_source_tables_are_frozen_with_exact_lookup_keys(self):
        self.assertEqual(self.data.source_revision, "9.0.5.0")
        self.assertEqual(len(self.data.alpha_parameters), 90)
        self.assertEqual(len(self.data.margules_interactions), 8)
        self.assertEqual(len(self.data.van_laar_interactions), 8)
        self.assertEqual(self.data.alpha("Methane").kappa1, -0.00159)
        self.assertEqual(self.data.alpha("methane").kappa1, 0.0)
        interaction = self.data.interaction("Acetone", "Cyclohexane", "margules")
        self.assertEqual((interaction.k12, interaction.k21), (0.0904, 0.1332))

    def test_prsv2_alpha_equation_matches_stated_source_vector(self):
        mixture = PRSV2Mixture(
            (self.catalog["Methane"], self.catalog["Ethane"]),
            (0.7, 0.3),
            self.data,
            "margules",
        )
        methane = self.catalog["Methane"]
        record = self.data.alpha("Methane")
        reduced_temperature = 180.0 / methane.critical_temperature.value
        root = math.sqrt(reduced_temperature)
        expected = (
            0.378893
            + 1.4897153 * methane.acentric_factor.value
            - 0.17131848 * methane.acentric_factor.value**2
            + 0.0196544 * methane.acentric_factor.value**3
            + (
                record.kappa1
                + record.kappa2
                * (record.kappa3 - reduced_temperature)
                * (1.0 - root)
            )
            * (1.0 + root)
            * (0.7 - reduced_temperature)
        )
        self.assertEqual(
            mixture.component_parameters(methane, 180.0).correction,
            expected,
        )

    def _assert_golden(
        self, stem: str, model_id: str, package_class: str, system_type: type
    ) -> None:
        golden = json.loads(
            (ROOT / f"tests/golden/{stem}.json").read_text(encoding="utf-8-sig")
        )
        repeat = json.loads(
            (ROOT / f"tests/golden/{stem}-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        self.assertEqual(golden, repeat)
        self.assertEqual(
            golden["source"]["property_package_class"],
            package_class,
        )
        inputs = golden["inputs"]
        system = create_thermo_system(
            model_id,
            compounds=tuple(self.catalog[name] for name in inputs["compounds"]),
            data=self.data,
        )
        self.assertIsInstance(system, system_type)
        for phase in ("liquid", "vapor"):
            state = system.state(
                tuple(inputs["composition"]),
                inputs["temperature_k"],
                inputs["pressure_pa"],
                phase,
            )
            for actual, expected in zip(
                state.fugacity_coefficients,
                golden["outputs"][f"{phase}_fugacity_coefficients"],
            ):
                self.assertTrue(math.isclose(actual, expected, rel_tol=5.0e-11))

    def test_prsv2_alpha_phase_states_match_repeatable_dwsim_golden(self):
        self._assert_golden(
            "prsv2-m-methane-ethane-state",
            PENG_ROBINSON_STRYJEK_VERA_2_MARGULES,
            "DWSIM.Thermodynamics.PropertyPackages.PRSV2PropertyPackage",
            PRSV2MargulesSystem,
        )

    def test_margules_phase_states_match_repeatable_dwsim_golden(self):
        self._assert_golden(
            "prsv2-m-acetone-cyclohexane-state",
            PENG_ROBINSON_STRYJEK_VERA_2_MARGULES,
            "DWSIM.Thermodynamics.PropertyPackages.PRSV2PropertyPackage",
            PRSV2MargulesSystem,
        )

    def test_van_laar_phase_states_match_repeatable_finite_dwsim_golden(self):
        self._assert_golden(
            "prsv2-vl-acetone-cyclohexane-state",
            PENG_ROBINSON_STRYJEK_VERA_2_VAN_LAAR,
            "DWSIM.Thermodynamics.PropertyPackages.PRSV2VLPropertyPackage",
            PRSV2VanLaarSystem,
        )

    def test_loader_rejects_duplicate_exact_alpha_key(self):
        document = json.loads(self.data_path.read_text(encoding="utf-8-sig"))
        document["alpha_parameters"][1]["compound"] = document["alpha_parameters"][0][
            "compound"
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaises(ValidationError):
                load_prsv2_data(path)


if __name__ == "__main__":
    unittest.main()
