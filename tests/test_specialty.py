import math
import sys
import unittest
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim.errors import ValidationError
from mesim.unitops.specialty import (
    DWSIM_GRAVITY_M_S2,
    DWSIM_WIND_POWER_COEFFICIENT,
    hydroelectric_turbine_power,
    solar_panel_power,
    wind_turbine_power,
)


ROOT = Path(__file__).parents[1]


def _saved_values(case_name, property_name):
    with ZipFile(ROOT / "tests" / case_name) as archive:
        xml_name = next(name for name in archive.namelist() if name.lower().endswith(".xml"))
        root = ElementTree.fromstring(archive.read(xml_name))
    return tuple(
        (element.text or "").strip()
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == property_name
    )


class SolarPanelParityTest(unittest.TestCase):
    def test_solar_panel_matches_fresh_dwsim_executable_save(self):
        irradiation = float(_saved_values(
            "u10-solar-panel.dwxmz", "ActualSolarIrradiation_kW_m2",
        )[0])
        area = float(_saved_values("u10-solar-panel.dwxmz", "PanelArea")[0])
        efficiency = float(_saved_values("u10-solar-panel.dwxmz", "PanelEfficiency")[0])
        count = int(_saved_values("u10-solar-panel.dwxmz", "NumberOfPanels")[0])
        reference = float(_saved_values("u10-solar-panel.dwxmz", "GeneratedPower")[0])

        result = solar_panel_power(irradiation, area, efficiency, count)

        self.assertTrue(math.isclose(result.generated_power_kw, reference, rel_tol=1e-14))
        self.assertIn("false", _saved_values("u10-solar-panel.dwxmz", "UseUserDefinedWeather"))

    def test_solar_panel_rejects_invalid_physical_inputs(self):
        with self.assertRaises(ValidationError):
            solar_panel_power(-1.0, 1.0, 15.0, 1)
        with self.assertRaises(ValidationError):
            solar_panel_power(1.0, 1.0, 101.0, 1)
        with self.assertRaises(ValidationError):
            solar_panel_power(1.0, 1.0, 15.0, 1.5)


class WindTurbineParityTest(unittest.TestCase):
    def test_wind_turbine_matches_official_dwsim_case(self):
        density = float(_saved_values("u10-wind-turbine.dwxmz", "AirDensity")[0])
        wind_speed = float(_saved_values("u10-wind-turbine.dwxmz", "WindSpeed_km_h")[0]) / 3.6
        area = float(_saved_values("u10-wind-turbine.dwxmz", "DiskArea")[0])
        efficiency = float(_saved_values("u10-wind-turbine.dwxmz", "Efficiency")[0])
        count = int(_saved_values("u10-wind-turbine.dwxmz", "NumberOfTurbines")[0])
        reference_maximum = float(_saved_values(
            "u10-wind-turbine.dwxmz", "MaximumTheoreticalPower",
        )[0])
        reference_generated = float(_saved_values(
            "u10-wind-turbine.dwxmz", "GeneratedPower",
        )[0])

        result = wind_turbine_power(density, wind_speed, area, efficiency, count)

        self.assertEqual(DWSIM_WIND_POWER_COEFFICIENT, 8.0 / 27.0)
        self.assertTrue(math.isclose(result.maximum_theoretical_power_kw, reference_maximum, rel_tol=1e-14))
        self.assertTrue(math.isclose(result.generated_power_kw, reference_generated, rel_tol=1e-14))
        self.assertTrue(math.isclose(result.rotor_diameter_m, math.sqrt(40.0 / math.pi)))

    def test_wind_turbine_rejects_invalid_physical_inputs(self):
        with self.assertRaises(ValidationError):
            wind_turbine_power(0.0, 10.0, 10.0, 80.0, 1)
        with self.assertRaises(ValidationError):
            wind_turbine_power(1.2, -1.0, 10.0, 80.0, 1)


class HydroelectricTurbineParityTest(unittest.TestCase):
    def test_hydroelectric_turbine_matches_official_dwsim_case_and_closes_energy(self):
        mass_flow = 997.048031971739
        density = 997.048498833022
        self.assertIn(str(mass_flow), _saved_values("u10-hydroelectric-turbine.dwxmz", "MassFlow"))
        self.assertIn(str(density), _saved_values("u10-hydroelectric-turbine.dwxmz", "density"))
        static_head = float(_saved_values("u10-hydroelectric-turbine.dwxmz", "StaticHead")[0])
        inlet_velocity = float(_saved_values("u10-hydroelectric-turbine.dwxmz", "InletVelocity")[0])
        outlet_velocity = float(_saved_values("u10-hydroelectric-turbine.dwxmz", "OutletVelocity")[0])
        efficiency = float(_saved_values("u10-hydroelectric-turbine.dwxmz", "Efficiency")[0])
        reference_velocity_head = float(_saved_values(
            "u10-hydroelectric-turbine.dwxmz", "VelocityHead",
        )[0])
        reference_total_head = float(_saved_values(
            "u10-hydroelectric-turbine.dwxmz", "TotalHead",
        )[0])
        reference_power = float(_saved_values(
            "u10-hydroelectric-turbine.dwxmz", "GeneratedPower",
        )[0])

        result = hydroelectric_turbine_power(
            mass_flow,
            density,
            static_head,
            inlet_velocity,
            outlet_velocity,
            efficiency,
        )

        self.assertEqual(DWSIM_GRAVITY_M_S2, 9.8)
        self.assertTrue(math.isclose(result.velocity_head_m, reference_velocity_head, rel_tol=1e-14))
        self.assertTrue(math.isclose(result.total_head_m, reference_total_head, rel_tol=1e-14))
        self.assertTrue(math.isclose(result.generated_power_kw, reference_power, rel_tol=1e-14))
        self.assertLess(abs(result.energy_balance_residual_kw), 1e-15)

    def test_hydroelectric_turbine_rejects_nonpositive_total_head(self):
        with self.assertRaises(ValidationError):
            hydroelectric_turbine_power(1000.0, 1000.0, -1.0, 0.0, 0.0, 75.0)


if __name__ == "__main__":
    unittest.main()
