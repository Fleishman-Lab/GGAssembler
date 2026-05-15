import os
import tempfile
import unittest
from typing import NamedTuple

from dawdlib.exact_construct import (
    Mutation,
    Variant,
    build_exact_variant_dnas,
    exact_variant_oligo_table,
    read_sequence_space_variants,
    variant_mutation_positions,
)
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.graph_maker import GraphMaker, make_no_degenerate_graph
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
        "K": "AAG",
        "L": "CTG",
        "T": "ACC",
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
        self.variants_csv = os.path.join(
            self.mthupo_dir, "MthUPO_150_variants.csv"
        )
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
        variant_dnas = build_exact_variant_dnas(
            wt_dna, variants, FakeCodonSelector()
        )

        self.assertEqual(variant_dnas["KLA"], "AAGCTTGCA")

    def test_variant_mutation_positions_uses_all_requested_mutations(self):
        variants = (
            Variant("v1", "v1", (Mutation(2, "A", "F"),)),
            Variant("v2", "v2", (Mutation(4, "A", "L"),)),
        )

        self.assertEqual(variant_mutation_positions(variants), [4, 5, 6, 10, 11, 12])

    def test_exact_variant_oligo_table_adds_only_variant_and_cassette_columns(self):
        path = [Gate(0, "AAAA"), Gate(8, "CCCC")]
        wt_dna = "A" * 30
        variant_dnas = {"v1": "A" * 30, "v2": "A" * 5 + "C" + "A" * 24}

        df = exact_variant_oligo_table(path, wt_dna, variant_dnas, "PRE", "SUF", "P")

        self.assertIn("variant_id", df.columns)
        self.assertIn("cassette_id", df.columns)
        self.assertEqual(df["cassette_id"].tolist(), [1, 1])
        self.assertEqual(df["const"].tolist(), [False, False])
        self.assertEqual(df["variant_id"].tolist(), ["v1", "v2"])


if __name__ == "__main__":
    unittest.main()
