import unittest

from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.graph_maker import (
    SOURCE,
    GraphMaker,
    create_default_weight_func,
    create_exact_variant_weight_func,
    make_exact_variant_graph,
    make_length_tier_cost,
    plan_exact_interval,
)
from dawdlib.golden_gate.utils import Requirements


class FakeGGData:
    def __init__(self, acceptable_overhangs):
        self.acceptable_overhangs = acceptable_overhangs

    def filter_self_binding_gates(self, _filter_gc):
        return self.acceptable_overhangs

    def overhangs_fidelity(self, _overhang1, _overhang2):
        return 1.0


def requirements(prefix="PRE", suffix="SUF"):
    return Requirements(
        min_oligo_length=0,
        max_oligo_length=100,
        min_const_oligo_length=0,
        min_efficiency=0.0,
        oligo_prefix=prefix,
        oligo_suffix=suffix,
        const_cost=0,
        min_fidelity=0.0,
        re_calc_lengths=False,
    )


class ExactVariantWeightTest(unittest.TestCase):
    def setUp(self):
        self.g1 = Gate(0, "AAAA")
        self.g2 = Gate(8, "CCCC")
        self.wt_dna = "A" * 30
        self.mutant_c = "A" * 5 + "C" + "A" * 24
        self.mutant_g = "A" * 5 + "G" + "A" * 24

    def test_default_individual_oligo_tiers(self):
        cost = make_length_tier_cost([(13, 60, 13), (61, 80, 17)])

        self.assertIsNone(cost(12))
        self.assertEqual(cost(13), 169)
        self.assertEqual(cost(60), 780)
        self.assertEqual(cost(61), 1037)
        self.assertEqual(cost(80), 1360)
        self.assertIsNone(cost(81))

    def test_length_tier_cost_uses_cheapest_match_and_explicit_fallback(self):
        cost = make_length_tier_cost(
            [(10, 20, 5), (15, 25, 2)], fallback=lambda length: length * 9
        )

        self.assertEqual(cost(16), 32)
        self.assertEqual(cost(30), 270)

    def test_length_cost_must_return_non_negative_integer_or_none(self):
        for invalid_cost in (-1, 1.5, "10", True):
            with self.subTest(invalid_cost=invalid_cost):
                with self.assertRaisesRegex(ValueError, "non-negative integer"):
                    plan_exact_interval(
                        self.g1,
                        self.g2,
                        self.wt_dna,
                        {"v1": self.mutant_c},
                        "PRE",
                        "SUF",
                        lambda _length, value=invalid_cost: value,
                    )

    def test_default_weight_func_empty_degenerate_table_keeps_const_cost(self):
        weight = create_default_weight_func({}, const_cost=7)

        self.assertEqual(weight(self.g1, self.g2), 7)

    def test_unique_mutant_payload_is_charged_once(self):
        plan = plan_exact_interval(
            self.g1,
            self.g2,
            self.wt_dna,
            {"v1": self.mutant_c, "v2": self.mutant_c},
            "PRE",
            "SUF",
            lambda length: length * 10,
        )

        self.assertEqual(len(plan.parts), 1)
        self.assertEqual(plan.parts[0].variant_ids, ("v1", "v2"))
        self.assertEqual(plan.parts[0].ordered_sequence, "PRE" + "AAAAACAAAAAA" + "SUF")
        self.assertEqual(plan.cost_cents, 180)

    def test_mixed_interval_sources_wt_and_two_mutants_independently(self):
        plan = plan_exact_interval(
            self.g1,
            self.g2,
            self.wt_dna,
            {"wt": self.wt_dna, "c": self.mutant_c, "g": self.mutant_g},
            "PRE",
            "SUF",
            lambda length: length * 10,
        )

        self.assertEqual(len(plan.parts), 3)
        self.assertEqual(
            [part.source_type for part in plan.parts],
            ["wt_pcr", "individual_oligo", "individual_oligo"],
        )
        self.assertEqual(plan.cost_cents, 520 + 180 + 180)

    def test_complete_ordered_length_controls_feasibility(self):
        cost = make_length_tier_cost([(12, 12, 1)])

        self.assertIsNone(
            plan_exact_interval(
                self.g1,
                self.g2,
                self.wt_dna,
                {"v1": self.mutant_c},
                "PRE",
                "SUF",
                cost,
            )
        )

    def test_individual_mode_has_no_retrieval_handles(self):
        plan = plan_exact_interval(
            self.g1,
            self.g2,
            self.wt_dna,
            {"v1": self.mutant_c},
            "PRE",
            "SUF",
            lambda length: length,
        )

        self.assertEqual(len(plan.parts[0].ordered_sequence), 18)
        self.assertEqual(plan.cost_cents, 18)

    def test_wt_direct_synthesis_is_optional(self):
        forced_pcr = plan_exact_interval(
            self.g1,
            self.g2,
            self.wt_dna,
            {"wt": self.wt_dna},
            "PRE",
            "SUF",
            lambda _length: 100,
            allow_wt_direct_synthesis=False,
        )
        cheapest = plan_exact_interval(
            self.g1,
            self.g2,
            self.wt_dna,
            {"wt": self.wt_dna},
            "PRE",
            "SUF",
            lambda _length: 100,
            allow_wt_direct_synthesis=True,
        )
        expensive_synthesis = plan_exact_interval(
            self.g1,
            self.g2,
            self.wt_dna,
            {"wt": self.wt_dna},
            "PRE",
            "SUF",
            lambda _length: 1000,
            allow_wt_direct_synthesis=True,
        )

        self.assertEqual(forced_pcr.parts[0].source_type, "wt_pcr")
        self.assertEqual(forced_pcr.cost_cents, 520)
        self.assertEqual(cheapest.parts[0].source_type, "individual_oligo")
        self.assertEqual(cheapest.cost_cents, 100)
        self.assertEqual(expensive_synthesis.parts[0].source_type, "wt_pcr")

    def test_synthesis_can_rescue_wt_part_too_short_for_pcr(self):
        plan = plan_exact_interval(
            self.g1,
            self.g2,
            self.wt_dna,
            {"wt": self.wt_dna},
            "PRE",
            "SUF",
            lambda _length: 200,
            min_const_oligo_length=20,
            allow_wt_direct_synthesis=True,
        )

        self.assertEqual(plan.parts[0].source_type, "individual_oligo")
        self.assertEqual(plan.cost_cents, 200)

    def test_terminal_wt_is_free_but_terminal_mutant_is_priced(self):
        plan = plan_exact_interval(
            SOURCE,
            self.g2,
            self.wt_dna,
            {"wt": self.wt_dna, "mutant": self.mutant_c},
            "PRE",
            "SUF",
            lambda length: length * 10,
        )

        self.assertEqual(plan.parts[0].source_type, "terminal_free")
        self.assertEqual(plan.parts[0].cost_cents, 0)
        self.assertEqual(plan.parts[1].ordered_sequence, "AAAAACAAAAAASUF")
        self.assertEqual(plan.parts[1].cost_cents, 150)

    def test_weight_function_returns_none_for_infeasible_edge(self):
        weight = create_exact_variant_weight_func(
            self.wt_dna,
            {"v1": self.mutant_c},
            prefix="PRE",
            suffix="SUF",
            segment_length_cost=lambda _length: None,
        )

        self.assertIsNone(weight(self.g1, self.g2))

    def test_exact_graph_uses_sequence_invariance_and_custom_filter(self):
        dna = "AAAACCCCGGGGTTTTAAAA"
        variant = dna[:5] + "A" + dna[6:]
        overhangs = {dna[i : i + 4] for i in range(len(dna) - 5)}
        graph, _source, _target = make_exact_variant_graph(
            GraphMaker(FakeGGData(overhangs)),
            dna.lower(),
            {"v1": variant.lower()},
            requirements(),
            segment_length_cost=lambda length: length,
            overhang_filter=lambda overhang: overhang != "AAAA",
        )

        regular_nodes = [node for node in graph if not node.src_or_target]
        self.assertNotIn("AAAA", {node.bps for node in regular_nodes})
        self.assertFalse(any(node.idx == 4 for node in regular_nodes))
        self.assertTrue(regular_nodes)


if __name__ == "__main__":
    unittest.main()
