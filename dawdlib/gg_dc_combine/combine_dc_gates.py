import json

import pandas as pd

from dawdlib.create_embl.embl_maker import df_to_gate_path
from dawdlib.gg_dc_combine.gg_dc_combine import (create_to_order_df,
                                                 parse_degenerate_codon_csv)
from dawdlib.golden_gate.utils import parse_dna


def combine_gate_path_deg_codons(
    dc_table_file: str,
    gate_path_file: str,
    dna_file: str,
    prefix: str,
    suffix: str,
    to_order_df_file: str,
):
    # dc_table_file = (
    #     "/home/labs/fleishman/jonathaw/for_others/200124_lihee/271/deg_table.csv"
    # )
    # gate_path_file = (
    #     "/home/labs/fleishman/jonathaw/for_others/200124_lihee/271/271_path_df.csv"
    # )
    # dna_file = (
    #     "/home/labs/fleishman/jonathaw/for_others/200124_lihee/271/271_2p48.fasta"
    # )
    # prefix = "CGTGCGGTCTCG"
    # suffix = "CGAGACCGCGCCGGGC"
    # to_order_df_file = (
    #     "/home/labs/fleishman/jonathaw/for_others/200124_lihee/271/to_order_temp.csv"
    # )

    dna = parse_dna(dna_file)
    dc_df = parse_degenerate_codon_csv(dc_table_file)
    path_df = pd.read_csv(gate_path_file)
    gate_path = df_to_gate_path(path_df)

    to_order_df = create_to_order_df(
        gate_path, deg_df=dc_df, dna=dna, prefix=prefix, suffix=suffix
    )
    to_order_df["wt"] = to_order_df["name"].str.contains("wt")
    to_order_df["gate1"] = to_order_df["gate1"].map(json.dumps)
    to_order_df["gate2"] = to_order_df["gate2"].map(json.dumps)
    to_order_df.to_csv(to_order_df_file)
