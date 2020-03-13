#!/usr/bin/env python3
import os
import unittest
from typing import List

import pandas as pd

from dawdlib.gg_dc_combine.gg_dc_combine import (
    create_dc_oligo,
    create_to_order_df,
    find_codons_for_oligo,
    find_oligos,
    parse_degenerate_codon_csv,
    parse_gg_segments_csv,
)
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.utils import find_dna_var_poss, parse_dna, parse_resfile


class GGDCCombneTest(unittest.TestCase):
    def setUp(self):
        self.here = os.path.dirname(__file__)
        self.gg_df = parse_gg_segments_csv(f"{self.here}/gg_segments.csv")
        self.dc_df = parse_degenerate_codon_csv(f"{self.here}/921.csv")
        self.dna = parse_dna(f"{self.here}/921_DNA_seq")

    def tearDown(self):
        pass

    def test_parse_gg_segments_csv(self):
        self.assertEqual(self.gg_df.iloc[0, 0], -1)
        self.assertEqual(self.gg_df.iloc[0, 1], "src")

        self.assertEqual(self.gg_df.iloc[1, 0], 79)
        self.assertEqual(self.gg_df.iloc[1, 1], "TTCG")

        self.assertEqual(self.gg_df.iloc[2, 0], 100)
        self.assertEqual(self.gg_df.iloc[2, 1], "ATGT")

    def test_parse_degenerate_codon_csv(self):
        self.assertEqual(self.dc_df.iloc[0, 0], 30)
        self.assertEqual(self.dc_df.iloc[0, 1], 88)
        self.assertEqual(self.dc_df.iloc[0, 2], ["A", "S", "T"])
        self.assertEqual(self.dc_df.iloc[0, 3], [1, 1, 1])
        self.assertEqual(self.dc_df.iloc[0, 4], "ASC")
        self.assertEqual(self.dc_df.iloc[0, 5], "GCA")

    def test_find_oligos(self):
        oligos = find_oligos(self.gg_df)
        self.assertListEqual(
            oligos,
            [
                (-1, 78),
                (79, 99),
                (100, 156),
                (157, 177),
                (178, 290),
                (291, 312),
                (313, 356),
            ],
        )

    def test_find_codons_for_oligo(self):
        oligos = find_oligos(self.gg_df)
        oligo_codons = find_codons_for_oligo(oligos[1], self.dc_df)
        self.assertEqual(
            oligo_codons,
            [
                [(88, "ASC"), (91, "AAA"), (97, "TWT")],
                [(88, "ASC"), (91, "ACC"), (97, "TWT")],
                [(88, "ASC"), (91, "CGT"), (97, "TWT")],
                [(88, "GCA"), (91, "AAA"), (97, "TWT")],
                [(88, "GCA"), (91, "ACC"), (97, "TWT")],
                [(88, "GCA"), (91, "CGT"), (97, "TWT")],
            ],
        )

    def test_create_dc_oligo(self):
        oligos = find_oligos(self.gg_df)
        oligo_codons = find_codons_for_oligo(oligos[1], self.dc_df)
        oligo_dna = create_dc_oligo(self.dna, oligo_codons[0], oligos[1])
        self.assertEqual(oligo_dna, "TCGACTTCASCAAATACTWTATGT")

    def test_create_to_order_df(self):
        deg_df = parse_degenerate_codon_csv(f"{self.here}/271_deg_table.csv")
        gg_df = parse_gg_segments_csv(f"{self.here}/271_path_df.csv")
        dna = parse_dna(f"{self.here}/271.fasta")
        prefix = "CGTGCGGTCTCG"
        suffix = "CGAGACCGCGCCGGGC"
        gg_list: List[Gate] = [
            Gate(idx=g.idx, bps=g.bps, req_primer=g.req_primer, syn_mut=g.syn_mut)
            for g in gg_df.itertuples()
        ]
        to_order_df = create_to_order_df(
            gate_path=gg_list, deg_df=deg_df, dna=dna, prefix=prefix, suffix=suffix
        )

        ref_df = pd.read_csv(f"{self.here}/271_to_order_df.csv")
        self.assertEqual(
            ref_df["oligo_dna"].tolist(), to_order_df["oligo_dna"].tolist()
        )
        print("AAAAAAAAaaa")


if __name__ == "__main__":
    unittest.main()
