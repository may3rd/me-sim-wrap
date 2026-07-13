import math
import sys
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import OutOfRangeError, ValidationError
from mesim.units import Quantity, convert, from_si, to_si


class UnitsTest(unittest.TestCase):
    def test_scale_conversions_cover_phase_two_dimensions(self):
        cases = (
            (1.0, "bar", "pressure", 100_000.0),
            (1.0, "kg/h", "mass_flow", 1 / 3600),
            (1.0, "kmol/h", "molar_flow", 1000 / 3600),
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
            to_si(1.0, "bar", "temperature")
        with self.assertRaises(ValidationError):
            to_si(float("nan"), "bar", "pressure")
        with self.assertRaises(ValidationError):
            to_si(True, "bar", "pressure")
        with self.assertRaises(OutOfRangeError):
            to_si(-274.0, "degC", "temperature")


if __name__ == "__main__":
    unittest.main()
