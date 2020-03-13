import filecmp
import os
import unittest

import pandas as pd


class IntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.here = os.path.dirname(__file__)
        # general input
        self.resfile_in = os.path.join(
            self.here, "resources", "271_2p48_reduced.resfile"
        )
        self.dna_in = os.path.join(self.here, "resources", "271_2p48.fasta")

        # partial inputs
        self.embl_in = os.path.join(self.here, "resources", "271_2p48.embl")
        self.embl_out = os.path.join(self.here, "test_output", "271_2p48.embl")
        self.embl_ref = os.path.join(self.here, "references", "w_features.gb")

        self.path_df_in = os.path.join(self.here, "resources", "271_path_df.csv")
        self.path_df_out = os.path.join(self.here, "test_output", "271_path_df.csv")
        self.path_df_ref = os.path.join(self.here, "references", "271_path_df.csv")

        self.deg_table_in = os.path.join(self.here, "resources", "deg_table.csv")
        self.deg_table_out = os.path.join(
            self.here, "test_output", "degenerate_table.csv"
        )
        self.deg_table_ref = os.path.join(
            self.here, "refernces", "degenerate_table.csv"
        )

        self.to_order_in = os.path.join(self.here, "resources", "to_order.csv")
        self.to_order_out = os.path.join(self.here, "test_output", "to_order.csv")
        self.to_order_ref = os.path.join(self.here, "references", "to_order.csv")

    def test_degenerate_table(self):
        from dawdlib.degenerate_dna.deg_table import generate_deg_csv

        generate_deg_csv(self.resfile_in, self.deg_table_out, "37762")
        self.assertTrue(filecmp.cmp(self.deg_table_out, self.deg_table_ref,))

    def test_create_goldengates(self):
        from dawdlib.golden_gate.find_gg import create_goldengates as cg_main

        cg_main(
            self.dna_in,
            self.deg_table_in,
            output_dir=os.path.join(self.here, "test_output"),
            min_var_oligo_length=20,
            max_var_oligo_length=79,
            min_const_oligo_length=20,
            no_solutions=1,
            min_oligos=None,
            max_oligos=None,
            gate_self_binding_min=2000,
            gate_crosstalk_max=1000,
        )

    def test_embl(self):
        from dawdlib.create_embl.create_embl import create_embl

        create_embl(
            deg_table_file=self.deg_table_in,
            embl_file=self.embl_in,
            path_file=self.path_df_in,
            out_embl_file=self.embl_out,
        )
        self.assertTrue(filecmp.cmp(self.embl_ref, self.embl_out,))

    def test_combine(self):
        from dawdlib.gg_dc_combine.combine_dc_gates import combine_gate_path_deg_codons

        combine_gate_path_deg_codons(
            dc_table_file=self.deg_table_in,
            gate_path_file=self.path_df_in,
            dna_file=self.dna_in,
            prefix="CGTGCGGTCTCG",
            suffix="CGAGACCGCGCCGGGC",
            to_order_df_file=self.to_order_out,
        )
        out_df = pd.read_csv(self.to_order_out)
        ref_df = pd.read_csv(self.to_order_ref)

        for col in out_df.columns:
            if "Unnamed" in col:
                continue
            self.assertTrue(all(out_df[col] == ref_df[col]))


if __name__ == "__main__":
    unittest.main()
