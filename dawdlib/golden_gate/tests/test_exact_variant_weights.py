import unittest

from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.graph_maker import (
    create_default_weight_func,
    create_exact_variant_weight_func,
    make_length_tier_cost,
)


class ExactVariantWeightTest(unittest.TestCase):
    def test_make_length_tier_cost_first_match_wins(self):
        cost = make_length_tier_cost([(10, 20, 2), (15, 25, 5)])

        self.assertEqual(cost(16), 32)
        self.assertEqual(cost(22), 110)
        self.assertEqual(cost(30), 30)

    def test_default_weight_func_empty_degenerate_table_keeps_const_cost(self):
        g1 = Gate(0, "AAAA")
        g2 = Gate(8, "CCCC")
        weight = create_default_weight_func({}, const_cost=7)

        self.assertEqual(weight(g1, g2), 7)

    def test_exact_weight_uses_const_cost_for_constant_edge(self):
        g1 = Gate(0, "AAAA")
        g2 = Gate(8, "CCCC")
        weight = create_exact_variant_weight_func(
            "A" * 30,
            {"v1": "A" * 30, "v2": "A" * 30},
            const_cost=11,
        )

        self.assertEqual(weight(g1, g2), 11)

    def test_exact_weight_deduplicates_unique_variant_segments(self):
        g1 = Gate(0, "AAAA")
        g2 = Gate(8, "CCCC")
        wt_dna = "A" * 30
        var_dna = "A" * 5 + "C" + "A" * 24
        weight = create_exact_variant_weight_func(
            wt_dna,
            {"v1": var_dna, "v2": var_dna},
            retrieval_handle_length=2,
            segment_length_cost=lambda length: length,
            oligo_addition=3,
            const_cost=11,
        )

        self.assertEqual(weight(g1, g2), 19)

    def test_exact_weight_includes_retrieval_handle_length(self):
        g1 = Gate(0, "AAAA")
        g2 = Gate(8, "CCCC")
        wt_dna = "A" * 30
        var_dna = "A" * 5 + "C" + "A" * 24
        weight = create_exact_variant_weight_func(
            wt_dna,
            {"v1": var_dna},
            retrieval_handle_length=5,
            segment_length_cost=lambda length: length * 2,
            oligo_addition=1,
            const_cost=11,
        )

        self.assertEqual(weight(g1, g2), 46)


if __name__ == "__main__":
    unittest.main()
