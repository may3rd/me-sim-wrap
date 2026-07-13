import math
import sys
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import OutOfRangeError, ValidationError
from mesim.units import Quantity, convert, from_si, to_si, unit_dimension


class UnitsTest(unittest.TestCase):
    def test_scale_conversions_cover_phase_two_dimensions(self):
        cases = (
            (1.0, "bar", "pressure", 100_000.0),
            (1.0, "kg/h", "mass_flow", 1 / 3600),
            (1.0, "kmol/h", "molar_flow", 1 / 3600),
            (1.0, "kJ", "energy", 1000.0),
            (1.0, "kW", "power", 1000.0),
            (1.0, "kJ/kg", "enthalpy", 1000.0),
            (1.0, "g/cm3", "density", 1000.0),
            (1.0, "cP", "viscosity", 0.001),
            (1.0, "W/m/K", "thermal_conductivity", 1.0),
            (1.0, "cm2", "area", 0.0001),
            (1.0, "L", "volume", 0.001),
            (1.0, "mm", "length", 0.001),
            (1.0, "ft/s", "velocity", 0.3048),
            (1.0, "kW/m2/K", "heat_transfer_coefficient", 1000.0),
            (100.0, "%", "dimensionless", 1.0),
        )

        for value, unit, dimension, expected in cases:
            with self.subTest(unit=unit):
                self.assertTrue(math.isclose(to_si(value, unit, dimension), expected, rel_tol=1e-12))

    def test_temperature_conversions_use_absolute_kelvin(self):
        self.assertEqual(to_si(273.15, "K", "temperature"), 273.15)
        self.assertTrue(math.isclose(to_si(0.0, "degC", "temperature"), 273.15, rel_tol=1e-12))
        self.assertTrue(math.isclose(to_si(32.0, "degF", "temperature"), 273.15, rel_tol=1e-12))
        self.assertTrue(math.isclose(from_si(273.15, "degC", "temperature"), 0.0, abs_tol=1e-12))

    def test_molar_enthalpy_uses_joules_per_kmol_as_si_base(self):
        cases = (
            ("J/kmol", 1.0),
            ("kJ/kmol", 1_000.0),
            ("J/mol", 1_000.0),
            ("kJ/mol", 1_000_000.0),
        )

        for unit, expected in cases:
            with self.subTest(unit=unit):
                self.assertEqual(to_si(1.0, unit, "molar_enthalpy"), expected)

    def test_process_unit_vectors_and_round_trips(self):
        bases = (
            ("delta_K", "temperature_difference"),
            ("s", "time"),
            ("kg", "mass"),
            ("kmol", "amount"),
            ("m3/s", "volumetric_flow"),
            ("J/kg/K", "heat_capacity"),
            ("J/kmol/K", "molar_heat_capacity"),
            ("J/kg/K", "entropy"),
            ("J/kmol/K", "molar_entropy"),
            ("m3/kg", "specific_volume"),
            ("m2/s", "kinematic_viscosity"),
            ("m2/s", "diffusivity"),
            ("N/m", "surface_tension"),
            ("Pa/m", "pressure_gradient"),
            ("kg/m2/s", "mass_flux"),
            ("W/m2", "heat_flux"),
            ("m/s2", "acceleration"),
            ("N", "force"),
        )
        cases = (
            (1.0, "delta_degF", "temperature_difference", 5.0 / 9.0),
            (0.0, "barg", "pressure", 101_325.0),
            (1.0, "h", "time", 3_600.0),
            (1.0, "lb", "mass", 0.45359237),
            (1.0, "mol", "amount", 0.001),
            (1.0, "m3/h", "volumetric_flow", 1.0 / 3_600.0),
            (1.0, "kJ/kg/K", "heat_capacity", 1_000.0),
            (1.0, "J/mol/K", "molar_heat_capacity", 1_000.0),
            (1.0, "kJ/kg/K", "entropy", 1_000.0),
            (1.0, "J/mol/K", "molar_entropy", 1_000.0),
            (1.0, "ft3/lb", "specific_volume", 0.028316846592 / 0.45359237),
            (1.0, "cSt", "kinematic_viscosity", 0.000001),
            (1.0, "cm2/s", "diffusivity", 0.0001),
            (1.0, "dyn/cm", "surface_tension", 0.001),
            (1.0, "bar/km", "pressure_gradient", 100.0),
            (1.0, "lb/ft2/s", "mass_flux", 0.45359237 / 0.09290304),
            (1.0, "Btu/h/ft2", "heat_flux", 1_055.05585262 / 3_600.0 / 0.09290304),
            (1.0, "ft/s2", "acceleration", 0.3048),
            (1.0, "lbf", "force", 4.4482216152605),
        )

        for unit, dimension in bases:
            with self.subTest(unit=unit):
                self.assertEqual(to_si(1.0, unit, dimension), 1.0)

        for value, unit, dimension, expected in cases:
            with self.subTest(unit=unit):
                actual = to_si(value, unit, dimension)
                self.assertTrue(math.isclose(actual, expected, rel_tol=1e-12))
                self.assertTrue(math.isclose(from_si(actual, unit, dimension), value, abs_tol=1e-12))

    def test_temperature_difference_is_not_absolute_temperature(self):
        self.assertEqual(to_si(20.0, "delta_degC", "temperature_difference"), 20.0)
        self.assertEqual(to_si(20.0, "degC", "temperature"), 293.15)
        self.assertEqual(unit_dimension("delta_K"), "temperature_difference")
        with self.assertRaises(ValidationError):
            to_si(20.0, "delta_degC", "temperature")

    def test_molar_flow_and_enthalpy_share_kmol_basis(self):
        flow_kmol_s = to_si(1.0, "kmol/s", "molar_flow")
        enthalpy_j_kmol = to_si(1.0, "kJ/kmol", "molar_enthalpy")

        self.assertEqual(flow_kmol_s * enthalpy_j_kmol, 1_000.0)

    def test_convert_round_trips_without_changing_requested_units(self):
        self.assertTrue(math.isclose(convert(1.0, "bar", "pressure", "psi"), 14.503773773, rel_tol=1e-9))
        self.assertTrue(math.isclose(convert(14.503773773, "psi", "pressure", "bar"), 1.0, rel_tol=1e-9))

    def test_quantity_preserves_input_and_is_immutable(self):
        quantity = Quantity.from_value(1.0, "bar", "pressure")

        self.assertEqual(quantity.value, 1.0)
        self.assertEqual(quantity.unit, "bar")
        self.assertEqual(quantity.si_value, 100_000.0)
        with self.assertRaises(FrozenInstanceError):
            quantity.value = 2.0

    def test_invalid_units_values_and_dimensions_are_rejected(self):
        with self.assertRaises(ValidationError):
            to_si(1.0, "not-a-unit", "pressure")
        with self.assertRaises(ValidationError):
            to_si(1.0, "gpm", "volumetric_flow")
        with self.assertRaises(ValidationError):
            to_si(1.0, "bar", "temperature")
        with self.assertRaises(ValidationError):
            to_si(float("nan"), "bar", "pressure")
        with self.assertRaises(ValidationError):
            to_si(True, "bar", "pressure")
        with self.assertRaises(ValidationError):
            to_si(1e308, "MPa", "pressure")
        with self.assertRaises(ValidationError):
            from_si(1e308, "cSt", "kinematic_viscosity")
        with self.assertRaises(OutOfRangeError):
            to_si(-274.0, "degC", "temperature")


if __name__ == "__main__":
    unittest.main()
