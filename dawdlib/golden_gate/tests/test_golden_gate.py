#!/usr/bin/env python3
import unittest
from dawdlib.golden_gate.graph_maker import GraphMaker, _is_node, _find_dna_var_poss
from dawdlib.golden_gate.gate_data import GGData


class GraphMakerTest(unittest.TestCase):
    def setUp(self):
        self.gg_data = GGData()
        self.graph_maker = GraphMaker(self.gg_data)

    def tearDown(self):
        pass

    def test_is_node(self):
        acceptable_fcws = self.gg_data.get_all_self_binding_gates()
        self.assertTrue(_is_node("AAAA", 10, acceptable_fcws, []))


if __name__ == '__main__':
    unittest.main()