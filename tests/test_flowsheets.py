import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.flowsheet import Connection, Flowsheet, UnitOperation, solve, topological_order
from mesim.streams import StreamState, flash_stream
from mesim.thermo.ideal import load_correlations
from mesim.unitops.basic import equilibrium_separator, heater, valve


ROOT = Path(__file__).parents[1]


class FlowsheetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compounds = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}
        cls.compounds = (compounds["Methane"], compounds["Ethane"])
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        cls.correlations = load_correlations(ROOT / "data/correlations/ideal-v1.json")

    @staticmethod
    def unit(tag, inputs=(), outputs=("out",), value=None):
        return UnitOperation(tag, inputs, outputs, lambda _: {"out": value} if not inputs else {"out": _["in"] + 1})

    def test_rejects_duplicate_tags_missing_ports_multiple_writers_and_cycles(self):
        feed = self.unit("feed", value=1)
        sink = self.unit("sink", ("in",))
        with self.assertRaises(ValidationError):
            topological_order(Flowsheet((feed, self.unit("feed", value=2)), ()))
        with self.assertRaises(ValidationError):
            topological_order(Flowsheet((feed, sink), ()))
        with self.assertRaises(ValidationError):
            topological_order(Flowsheet((feed, sink), (Connection("feed", "missing", "sink", "in"),)))
        other = self.unit("other", value=2)
        with self.assertRaises(ValidationError):
            topological_order(Flowsheet((feed, other, sink), (
                Connection("feed", "out", "sink", "in"), Connection("other", "out", "sink", "in"),
            )))
        first = self.unit("first", ("in",))
        second = self.unit("second", ("in",))
        with self.assertRaises(ValidationError):
            topological_order(Flowsheet((first, second), (
                Connection("first", "out", "second", "in"), Connection("second", "out", "first", "in"),
            )))

    def test_orders_and_solves_acyclic_graph_without_mutating_prior_results(self):
        feed = self.unit("feed", value=1)
        add_a = self.unit("add-a", ("in",))
        add_b = self.unit("add-b", ("in",))
        flowsheet = Flowsheet((add_b, feed, add_a), (
            Connection("feed", "out", "add-a", "in"), Connection("add-a", "out", "add-b", "in"),
        ))

        self.assertEqual(topological_order(flowsheet), ("feed", "add-a", "add-b"))
        result = solve(flowsheet)
        self.assertEqual(result.values[("add-b", "out")], 3)
        with self.assertRaises(TypeError):
            result.values[("add-b", "out")] = 4

    def test_executes_feed_heater_valve_separator_without_partial_results(self):
        def feed(_):
            return {"out": flash_stream(
                StreamState(180.0, 500_000.0, 2.0, ("Methane", "Ethane"), (0.7, 0.3)),
                self.compounds, self.interactions, self.correlations,
            )}

        flowsheet = Flowsheet(
            (
                UnitOperation("separator", ("in",), ("liquid", "vapor"), lambda values: {
                    "liquid": equilibrium_separator(values["in"], self.compounds, self.correlations).liquid,
                    "vapor": equilibrium_separator(values["in"], self.compounds, self.correlations).vapor,
                }),
                UnitOperation("valve", ("in",), ("out",), lambda values: {"out": valve(
                    values["in"], self.compounds, self.interactions, self.correlations, 300_000.0, (140.0, 240.0),
                )}),
                UnitOperation("feed", (), ("out",), feed),
                UnitOperation("heater", ("in",), ("out",), lambda values: {"out": heater(
                    values["in"], self.compounds, self.interactions, self.correlations, 190.0,
                ).outlet}),
            ),
            (
                Connection("feed", "out", "heater", "in"),
                Connection("heater", "out", "valve", "in"),
                Connection("valve", "out", "separator", "in"),
            ),
        )

        result = solve(flowsheet)
        liquid, vapor = result.values[("separator", "liquid")], result.values[("separator", "vapor")]
        self.assertIsNotNone(liquid)
        self.assertIsNotNone(vapor)
        self.assertEqual(liquid.molar_flow_kmol_s + vapor.molar_flow_kmol_s, 2.0)
        self.assertEqual(liquid.pressure_pa, vapor.pressure_pa)


if __name__ == "__main__":
    unittest.main()
