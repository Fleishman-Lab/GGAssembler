#!/usr/bin/env python3
import os
import unittest
from dawdlib.golden_gate.graph_maker import GraphMaker, _is_node
from dawdlib.golden_gate.utils import find_dna_var_poss
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.utils import parse_dna


class GraphMakerTest(unittest.TestCase):
    def setUp(self):
        self.gg_data = GGData()
        self.graph_maker = GraphMaker(self.gg_data)
        self.dna = parse_dna(f"{os.path.dirname(__file__)}/input_921/921_DNA_seq")

    def tearDown(self):
        pass

    def test_is_node(self):
        acceptable_fcws = self.gg_data.get_all_self_binding_gates()
        self.assertTrue(not _is_node("AAAA", 10, acceptable_fcws, []))
        self.assertTrue(_is_node("AAAC", 10, acceptable_fcws, list(range(100))))


if __name__ == '__main__':
    unittest.main()