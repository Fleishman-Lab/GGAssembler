from typing import List

import pandas as pd
from Bio import SeqFeature, SeqIO

from dawdlib.degenerate_dna.utils import parse_degenerate_codon_csv
from dawdlib.embl_utils.embl_maker import create_dc_features, create_path_features
from dawdlib.golden_gate.utils import gate_df_list


def create_embl(
    deg_table_file: str, embl_file: str, path_file: str, out_embl_file: str
):
    """

    Args:
        deg_table_file (str):
        embl_file (str):
        path_file (str):
        out_embl_file (str):
    """
    seq_record = SeqIO.read(embl_file, "embl")

    features: List[SeqFeature.SeqFeature] = []

    if deg_table_file:
        deg_parsed_df = parse_degenerate_codon_csv(deg_table_file)
        features.extend(create_dc_features(deg_parsed_df))

    if path_file:
        path_df = pd.read_csv(path_file)
        gate_path = gate_df_list(path_df)
        features.extend(create_path_features(gate_path))

    seq_record.features.extend(features)
    with open(out_embl_file, "w+") as fout:
        SeqIO.write(seq_record, fout, "gb")
