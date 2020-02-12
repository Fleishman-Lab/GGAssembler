from typing import List

import pandas as pd
from Bio import SeqFeature, SeqIO

from dawdlib.create_embl.embl_maker import (create_dc_features, create_path_features, df_to_gate_path,
                                            parse_degenerate_codon_csv)
from dawdlib.golden_gate.utils import parse_dna


def main():
    deg_table_file = "/home/labs/fleishman/jonathaw/for_others/200124_lihee/271/deg_table.csv"
    embl_file = "/home/labs/fleishman/jonathaw/for_others/200124_lihee/271/271_2p48.embl"
    dna_file = ""
    path_file = "/home/labs/fleishman/jonathaw/for_others/200124_lihee/271/271_path_df.csv"
    out_embl = "/home/labs/fleishman/jonathaw/for_others/200124_lihee/271/out.embl"

    seq_record = SeqIO.read(embl_file, "embl")
    # dna = parse_dna(dna_file)

    features: List[SeqFeature.SeqFeature] = []

    if deg_table_file:
        deg_parsed_df = parse_degenerate_codon_csv(deg_table_file)
        features.extend(create_dc_features(deg_parsed_df))

    if path_file:
        path_df = pd.read_csv(path_file)
        gate_path = df_to_gate_path(path_df)
        features.extend(create_path_features(gate_path))

    seq_record.features.extend(features)
    with open(out_embl, "w+") as fout:
        SeqIO.write(seq_record, fout, "gb")


if __name__ == "__main__":
    main()