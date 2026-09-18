import os
import tempfile
import unittest
from typing import NamedTuple

import networkx as nx

from dawdlib.exact_construct import (
    Mutation,
    Variant,
    build_exact_variant_dnas,
    exact_variant_tables,
    read_sequence_space_variants,
    variant_mutation_positions,
)
from dawdlib.dijkstra.len_limit import all_shortest_paths
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gate_restriction import create_path_validator
from dawdlib.golden_gate.graph_maker import (
    SOURCE,
    TARGET,
    GraphMaker,
    make_exact_variant_graph,
    make_length_tier_cost,
    make_no_degenerate_graph,
    plan_exact_interval,
)
from dawdlib.golden_gate.utils import RequirementsFactory, parse_dna


class FakeCodon(NamedTuple):
    codon: str
    probability: float
    cai: float


class FakeAminoAcid(NamedTuple):
    amino_acid: str
    codons: tuple


class FakeSelection(NamedTuple):
    amino_acids: tuple
    ambiguous_codons: tuple


class FakeCodonSelector:
    codons = {
        "A": "GCG",
        "F": "TTC",
        "G": "GGC",
        "I": "ATC",
        "K": "AAG",
        "L": "CTG",
        "M": "ATG",
        "S": "AGC",
        "T": "ACC",
        "Y": "TAC",
    }

    def optimise_codons(self, amino_acids):
        aa = amino_acids[0]
        codon = self.codons[aa]
        return FakeSelection(
            amino_acids=(FakeAminoAcid(aa, (FakeCodon(codon, 1.0, 1.0),)),),
            ambiguous_codons=(codon,),
        )


class ExactConstructTest(unittest.TestCase):
    def setUp(self):
        self.mthupo_dir = os.path.join(
            os.path.dirname(__file__), "..", "..", "example", "MthUPO"
        )
        self.variants_csv = os.path.join(self.mthupo_dir, "MthUPO_150_variants.csv")
        self.resfile = os.path.join(self.mthupo_dir, "MthUPO_Sequence_Space.resfile")
        self.wt_dna = parse_dna(os.path.join(self.mthupo_dir, "MthUPO_wt.txt"))

    def test_no_degenerate_graph_builds(self):
        reqs = RequirementsFactory(
            4,
            40,
            4,
            0.1,
            "",
            "",
            1,
            min_fidelity=0.0,
        )
        graph, source, target = make_no_degenerate_graph(
            GraphMaker(GGData()), "ATG" * 20, reqs
        )

        self.assertIn(source, graph.nodes)
        self.assertIn(target, graph.nodes)
        self.assertTrue(graph.edges)

    def test_read_sequence_space_variants_parses_mthupo_inputs(self):
        variants = read_sequence_space_variants(
            self.variants_csv, self.resfile, self.wt_dna
        )

        self.assertEqual(len(variants), 150)
        self.assertEqual(variants[0].variant_id, "KLLLLLLATF")
        self.assertEqual(variants[0].total_score, -776.16)

    def test_sequence_space_variant_rejects_wrong_length_sequence_name(self):
        with tempfile.NamedTemporaryFile("w", suffix=".csv") as handle:
            handle.write("sequence_name,total_score\nKLL,-1\n")
            handle.flush()

            with self.assertRaisesRegex(ValueError, "expected 10"):
                read_sequence_space_variants(handle.name, self.resfile, self.wt_dna)

    def test_sequence_space_variant_rejects_disallowed_amino_acid(self):
        with tempfile.NamedTemporaryFile("w", suffix=".csv") as handle:
            handle.write("sequence_name,total_score\nZLLLLLLATF,-1\n")
            handle.flush()

            with self.assertRaisesRegex(ValueError, "not allowed"):
                read_sequence_space_variants(handle.name, self.resfile, self.wt_dna)

    def test_first_mthupo_variant_maps_resfile_positions_in_order(self):
        variant = read_sequence_space_variants(
            self.variants_csv, self.resfile, self.wt_dna
        )[0]

        self.assertEqual(
            variant.mutations,
            (
                Mutation(59, "F", "K"),
                Mutation(63, "F", "L"),
                Mutation(153, "A", "L"),
                Mutation(154, "F", "L"),
                Mutation(156, "Y", "L"),
                Mutation(157, "G", "A"),
                Mutation(159, "S", "T"),
                Mutation(161, "A", "F"),
            ),
        )

    def test_sequence_space_variant_uses_translated_wt_dna_for_from_aa(self):
        variant = read_sequence_space_variants(
            self.variants_csv, self.resfile, self.wt_dna
        )[0]

        self.assertEqual(variant.mutations[0], Mutation(59, "F", "K"))

    def test_exact_variant_dna_generation_from_mutation_records(self):
        wt_dna = "TTTCTTGCA"
        variants = (Variant("KLA", "KLA", (Mutation(1, "F", "K"),)),)
        variant_dnas = build_exact_variant_dnas(wt_dna, variants, FakeCodonSelector())

        self.assertEqual(variant_dnas["KLA"], "AAGCTTGCA")

    def test_exact_variant_dna_generation_rejects_duplicate_ids(self):
        variants = (
            Variant("same", "KLA", (Mutation(1, "F", "K"),)),
            Variant("same", "FLA", ()),
        )

        with self.assertRaisesRegex(ValueError, "IDs must be unique"):
            build_exact_variant_dnas("TTTCTTGCA", variants, FakeCodonSelector())

    def test_variant_mutation_positions_uses_all_requested_mutations(self):
        variants = (
            Variant("v1", "v1", (Mutation(2, "A", "F"),)),
            Variant("v2", "v2", (Mutation(4, "A", "L"),)),
        )

        self.assertEqual(variant_mutation_positions(variants), [4, 5, 6, 10, 11, 12])

    def test_exact_variant_tables_deduplicate_parts_and_report_summary(self):
        wt_dna = "A" * 30
        variant_dnas = {"v1": "A" * 30, "v2": "A" * 5 + "C" + "A" * 24}
        gate1 = Gate(0, "AAAA")
        gate2 = Gate(8, "CCCC")
        target = TARGET._replace(idx=len(wt_dna) + 1)
        path = [SOURCE, gate1, gate2, target]
        graph = nx.DiGraph()
        for first, second in zip(path[:-1], path[1:]):
            plan = plan_exact_interval(
                first,
                second,
                wt_dna,
                variant_dnas,
                "PRE",
                "SUF",
                lambda length: length * 10,
            )
            graph.add_edge(first, second, weight=plan.cost_cents, exact_plan=plan)

        output = exact_variant_tables(
            path, graph, variant_dnas, fidelity_metrics={"NEB": 0.95}
        )

        self.assertEqual(len(output.parts), 4)
        self.assertEqual(len(output.bom), 6)
        internal_bom = output.bom[output.bom.interval_id == "I2"]
        self.assertEqual(internal_bom.part_id.nunique(), 2)
        self.assertEqual(output.summary["total_cost_cents"], 700)
        self.assertEqual(output.summary["unique_synthesized_non_wt_oligos"], 1)
        self.assertEqual(output.summary["wt_pcr_fragments"], 1)
        self.assertEqual(output.summary["terminal_free_parts"], 2)
        self.assertEqual(output.summary["total_ordered_oligo_nt"], 18)
        self.assertEqual(output.summary["total_primer_nt"], 40)
        self.assertEqual(output.summary["total_fragment_usages"], 6)
        self.assertEqual(output.summary["average_non_wt_parts_per_variant"], 0.5)
        self.assertEqual(output.summary["average_non_common_parts_per_variant"], 1.0)
        self.assertEqual(output.summary["golden_gate_fragments_per_construct"], 3)
        self.assertEqual(output.summary["golden_gate_junctions"], 2)
        self.assertEqual(output.summary["fidelity"], {"NEB": 0.95})

    def test_existing_path_search_selects_lower_acquisition_cost(self):
        graph = nx.DiGraph()
        graph.add_weighted_edges_from(
            [
                ("source", "expensive", 50),
                ("expensive", "target", 50),
                ("source", "cheap1", 10),
                ("cheap1", "cheap2", 10),
                ("cheap2", "target", 10),
            ]
        )

        paths = list(
            all_shortest_paths(graph, "source", "target", len_cutoff=4, weight="weight")
        )

        self.assertEqual(
            min(paths, key=lambda path_and_cost: path_and_cost[1]),
            (["source", "cheap1", "cheap2", "target"], 30),
        )

    def test_mthupo_exact_graph_path_and_outputs(self):
        variants = read_sequence_space_variants(
            self.variants_csv, self.resfile, self.wt_dna
        )
        variant_dnas = build_exact_variant_dnas(
            self.wt_dna, variants, FakeCodonSelector(), enzymes=["BsaI"]
        )
        reqs = RequirementsFactory(
            min_oligo_length=0,
            max_oligo_length=10_000,
            min_const_oligo_length=15,
            min_efficiency=0.25,
            oligo_prefix="GACATTGGTCTCA",
            oligo_suffix="TGAGACCAACGACGCCGTACTCTTTGTCAAC",
            const_cost=0,
            min_fidelity=0.1,
            filter_gc_overhangs=False,
        )
        ggdata = GGData(
            temperature=reqs.gg_temp,
            hours=reqs.gg_hours,
            min_efficiency=reqs.min_efficiency,
            min_fidelity=reqs.min_fidelity,
        )
        graph, source, target = make_exact_variant_graph(
            GraphMaker(ggdata),
            self.wt_dna,
            variant_dnas,
            reqs,
            segment_length_cost=make_length_tier_cost([(13, 60, 13), (61, 80, 17)]),
        )
        junctions = (
            (169, "CCAG"),
            (180, "TTTG"),
            (191, "ATTA"),
            (248, "CAAC"),
            (259, "AACA"),
            (452, "AGAT"),
            (483, "TATG"),
        )
        nodes = {(node.idx, node.bps): node for node in graph if not node.src_or_target}
        path = [source, *(nodes[junction] for junction in junctions), target]

        self.assertTrue(create_path_validator(ggdata)(path))
        cost_cents = sum(
            graph.edges[first, second]["weight"]
            for first, second in zip(path[:-1], path[1:])
        )
        output = exact_variant_tables(path, graph, variant_dnas)

        self.assertEqual(cost_cents, 200_785)
        self.assertEqual(output.summary["total_cost_cents"], cost_cents)
        self.assertEqual(output.summary["unique_physical_parts"], 158)
        self.assertEqual(len(output.bom), 1_200)


if __name__ == "__main__":
    unittest.main()
