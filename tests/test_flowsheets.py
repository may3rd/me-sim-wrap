import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.flowsheet import Connection, Flowsheet, UnitOperation, solve, topological_order


class FlowsheetTest(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
