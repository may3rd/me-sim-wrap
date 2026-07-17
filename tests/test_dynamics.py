import csv
import math
import sys
import unittest
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim.errors import ValidationError
from mesim.unitops.dynamics import (
    AlgebraicConstraint,
    HeatExchangerState,
    HoldupRates,
    HoldupState,
    PIDConfig,
    StateEvent,
    TankLevelControlConfig,
    adaptive_explicit_ode,
    advance_holdup,
    dynamic_cstr_step,
    fixed_step_explicit_euler,
    lumped_heat_exchanger_step,
    simulate_dwsim_tank_level_control,
    validate_algebraic_constraints,
)


ROOT = Path(__file__).parents[1]
DWSIM_LEVEL_ABS_TOL_M = 1.0e-4
DWSIM_VALVE_ABS_TOL_PERCENT = 1.3e-2


def _local_name(element):
    return element.tag.rsplit("}", 1)[-1]


def _saved_values(root, name):
    return [(element.text or "").strip() for element in root.iter() if _local_name(element) == name]


class DynamicFoundationTest(unittest.TestCase):
    def test_holdup_closes_component_and_energy_accumulation(self):
        result = advance_holdup(
            HoldupState((2.0, 1.0), 2_000_000.0),
            HoldupRates(
                component_in_kmol_s=(0.2, 0.1),
                component_out_kmol_s=(0.1, 0.08),
                component_generation_kmol_s=(-0.01, 0.01),
                energy_in_w=50_000.0,
                energy_out_w=40_000.0,
                heat_input_w=-2_000.0,
                energy_generation_w=1_000.0,
            ),
            5.0,
        )

        self.assertEqual(result.state.component_amounts_kmol, (2.45, 1.15))
        self.assertEqual(result.state.internal_energy_j, 2_045_000.0)
        self.assertLess(max(map(abs, result.component_balance_residuals_kmol)), 5e-16)
        self.assertEqual(result.energy_balance_residual_j, 0.0)

    def test_algebraic_initialization_is_explicit_and_scaled(self):
        records = validate_algebraic_constraints(
            (AlgebraicConstraint("volume closure", 1e-8, 2.0),), 1e-7,
        )
        self.assertEqual(records[0].name, "volume closure")
        with self.assertRaisesRegex(ValidationError, "inconsistent algebraic initialization"):
            validate_algebraic_constraints(
                (AlgebraicConstraint("pressure closure", 1e-3, 1.0),), 1e-6,
            )

    def test_fixed_step_mode_is_reproducible_and_applies_aligned_events(self):
        arguments = (
            lambda _time, state: (-state[0],),
            (1.0,),
            1.0,
            0.25,
            (StateEvent(0.5, 0, 2.0),),
        )
        first = fixed_step_explicit_euler(*arguments)
        second = fixed_step_explicit_euler(*arguments)

        self.assertEqual(first, second)
        self.assertEqual(first.times_s, (0.0, 0.25, 0.5, 0.75, 1.0))
        self.assertEqual(first.states[2], (2.0,))
        self.assertEqual(first.states[-1], (1.125,))
        with self.assertRaisesRegex(ValidationError, "align"):
            fixed_step_explicit_euler(
                lambda _time, state: state, (1.0,), 1.0, 0.25, (StateEvent(0.3, 0, 2.0),),
            )

    def test_adaptive_mode_integrates_reduced_ode_and_rejects_dae_claims(self):
        trajectory = adaptive_explicit_ode(
            lambda _time, state: (-state[0],),
            (1.0,),
            1.0,
            output_times_s=(0.0, 0.5, 1.0),
            relative_tolerance=1e-10,
            absolute_tolerance=1e-12,
        )
        self.assertTrue(math.isclose(trajectory.states[-1][0], math.exp(-1.0), rel_tol=2e-10))
        with self.assertRaisesRegex(ValidationError, "IDA-capable"):
            adaptive_explicit_ode(
                lambda _time, state: state,
                (1.0,),
                1.0,
                algebraic_equations_present=True,
            )

    def test_lumped_heat_exchanger_closes_total_energy_each_step(self):
        result = lumped_heat_exchanger_step(
            HeatExchangerState(400.0, 300.0),
            hot_inlet_temperature_k=420.0,
            cold_inlet_temperature_k=290.0,
            hot_mass_flow_kg_s=2.0,
            cold_mass_flow_kg_s=3.0,
            hot_holdup_kg=100.0,
            cold_holdup_kg=150.0,
            hot_heat_capacity_j_kg_k=2_000.0,
            cold_heat_capacity_j_kg_k=4_000.0,
            ua_w_k=500.0,
            step_s=2.0,
        )

        self.assertEqual(result.transferred_heat_w, 50_000.0)
        self.assertEqual(result.state.hot_temperature_k, 400.3)
        self.assertTrue(math.isclose(result.state.cold_temperature_k, 299.76666666666665))
        self.assertEqual(result.total_energy_balance_residual_j, 0.0)

    def test_dynamic_cstr_uses_one_kmol_and_j_basis(self):
        result = dynamic_cstr_step(
            HoldupState((5.0, 0.0), 1_000_000.0),
            component_in_kmol_s=(0.2, 0.0),
            component_out_kmol_s=(0.1, 0.0),
            stoichiometry=((-1.0,), (1.0,)),
            reaction_extents_kmol_s=(0.02,),
            energy_in_w=20_000.0,
            energy_out_w=15_000.0,
            reaction_enthalpies_j_kmol=(-100_000.0,),
            heat_input_w=-1_000.0,
            step_s=10.0,
        )

        self.assertTrue(math.isclose(result.state.component_amounts_kmol[0], 5.8))
        self.assertTrue(math.isclose(result.state.component_amounts_kmol[1], 0.2))
        self.assertEqual(result.state.internal_energy_j, 1_060_000.0)
        self.assertLess(max(map(abs, result.component_balance_residuals_kmol)), 1e-15)
        self.assertEqual(result.energy_balance_residual_j, 0.0)


class DwsimTankLevelParityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.primary_path = ROOT / "tests/golden/u9-water-tank-level-control-dwsim.txt"
        cls.repeat_path = ROOT / "tests/golden/u9-water-tank-level-control-dwsim-repeat.txt"
        with cls.primary_path.open(encoding="utf-8-sig", newline="") as source:
            cls.reference = tuple(csv.DictReader(source, delimiter="\t"))

    def test_dwsim_capture_is_repeatable_and_has_no_missing_records(self):
        self.assertEqual(
            self.primary_path.read_bytes(),
            self.repeat_path.read_bytes(),
        )
        self.assertEqual(len(self.reference), 121)
        self.assertEqual(tuple(self.reference[0]), ("Time (s)", "Level (m)", "V1 (%)", "V2 (%)"))
        self.assertEqual(float(self.reference[0]["Time (s)"]), 0.0)
        self.assertEqual(float(self.reference[-1]["Time (s)"]), 600.0)
        for row in self.reference:
            for value in row.values():
                self.assertTrue(math.isfinite(float(value)))

    def test_saved_dwsim_case_freezes_schedule_vessel_valve_and_pid_inputs(self):
        path = ROOT / "tests/u9-water-tank-level-control.dwxmz"
        with ZipFile(path) as archive:
            xml_name = next(name for name in archive.namelist() if name.lower().endswith(".xml"))
            root = ElementTree.fromstring(archive.read(xml_name))

        self.assertIn("5000", _saved_values(root, "IntegrationStep"))
        self.assertIn("604000", _saved_values(root, "Duration"))
        self.assertIn("false", _saved_values(root, "UseCurrentStateAsInitial"))
        self.assertIn("2", _saved_values(root, "Volume"))
        self.assertIn("400", _saved_values(root, "Kv"))
        self.assertIn("33.8440722851719", _saved_values(root, "OpeningPct"))
        self.assertIn("119.455297118019", _saved_values(root, "Kp"))
        self.assertIn("4.52783924040604", _saved_values(root, "Ki"))
        self.assertIn("16.3382185733144", _saved_values(root, "Kd"))
        self.assertIn("1.7", _saved_values(root, "SPValue"))

    def test_python_tank_pid_trajectory_matches_every_dwsim_record(self):
        density = 997.060396254973
        inlet_volume_flow = 10.0 / density
        config = TankLevelControlConfig(
            duration_s=600.0,
            step_s=5.0,
            tank_volume_m3=2.0,
            tank_height_m=2.0,
            liquid_density_kg_m3=density,
            inlet_mass_flow_kg_s=10.0,
            initial_contents_volume_m3=inlet_volume_flow,
            initial_outlet_opening_percent=33.8440722851719,
            inlet_opening_percent=50.0,
            tank_base_pressure_pa=101325.0,
            downstream_pressure_pa=109909.260185509,
            gravity_m_s2=9.8,
            valve_kv=400.0,
            pid=PIDConfig(
                proportional_gain=119.455297118019,
                integral_gain=4.52783924040604,
                derivative_gain=16.3382185733144,
                setpoint=1.7,
                output_scale=1.7,
                integral_guard=20.0,
                reverse_acting=True,
                output_min=-1000.0,
                output_max=1000.0,
            ),
        )

        trajectory = simulate_dwsim_tank_level_control(config)

        self.assertEqual(len(trajectory.points), len(self.reference))
        self.assertLess(max(map(abs, trajectory.mass_balance_residuals_kg)), 4e-13)
        for point, row in zip(trajectory.points, self.reference, strict=True):
            self.assertEqual(point.time_s, float(row["Time (s)"]))
            self.assertEqual(point.inlet_opening_percent, float(row["V1 (%)"]))
            self.assertLessEqual(
                abs(point.liquid_level_m - float(row["Level (m)"])),
                DWSIM_LEVEL_ABS_TOL_M,
            )
            self.assertLessEqual(
                abs(point.outlet_opening_percent - float(row["V2 (%)"])),
                DWSIM_VALVE_ABS_TOL_PERCENT,
            )


if __name__ == "__main__":
    unittest.main()
