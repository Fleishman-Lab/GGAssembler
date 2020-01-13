from itertools import product
from typing import List, Tuple

import numpy as np
import pandas as pd


def parse_gg_segments_csv(csv_file: str) -> pd.DataFrame:
    return pd.read_csv(csv_file, index_col=False)


def parse_degenerate_codon_csv(csv_file: str) -> pd.DataFrame:
    return pd.read_csv(
        csv_file,
        index_col=False,
        na_values="NaN",
        converters={
            "ENCODED_AAS": lambda x: x.strip("[]").replace("'", "").split(", "),
            "ENCODED_COUNT": lambda x: [
                int(a) for a in x.strip("[]").replace("'", "").split(", ")
            ],
        },
    )


def find_oligos(gg_df: pd.DataFrame) -> List[Tuple[int, int]]:
    oligos: List[Tuple[int, int]] = []
    for ind1, ind2 in zip(gg_df.index[:-1], gg_df.index[1:]):
        oligos.append((gg_df.loc[ind1, "idx"], gg_df.loc[ind2, "idx"] - 1))
    return oligos


def find_codons_for_oligo(
    oligo: Tuple[int, int], dc_df: pd.DataFrame
) -> List[List[Tuple[int, str]]]:

    oligo_codons: List[List[Tuple[int, str]]] = []
    sub_df = dc_df.loc[
        dc_df["DNA_POS"].between(oligo[0], oligo[1]),
        ["DNA_POS"] + [c for c in dc_df.columns if c.startswith("AMBIGUOUS_CODONS")],
    ]
    pos_codon: List[List[Tuple[int, str]]] = []
    for _, row in sub_df.iterrows():
        pos_codon.append([])
        for col in sub_df.columns[1:]:
            if row[col] is np.nan:
                continue
            pos_codon[-1].append((row["DNA_POS"], row[col]))

    for prod in product(*pos_codon):
        combination: List[Tuple[int, str]] = []
        for (pos, codon) in prod:
            combination.append((pos, codon))
        oligo_codons.append(combination)
    return oligo_codons

def create_dc_oligo(dna: str, pos_codons: List[Tuple[int, str]], oligo: Tuple[int, int]) -> str:
    dna_copy = dna
    for (pos, codon) in pos_codons:
        dna_copy = dna_copy[:pos] + codon + dna_copy[pos+4:]
    return dna_copy[oligo[0]: oligo[1]+3]


def create_all_dc_oligos(dna: str, gg_df: pd.DataFrame, dc_df: pd.DataFrame):
    oligos = find_oligos(gg_df)
    for oligo in oligos:
        oligo_codons = find_codons_for_oligo(oligo, dc_df)
        for oligo_codon in oligo_codons:
            oligo_dna = create_dc_oligo(dna, oligo_codon, oligo)
            print(oligo, oligo_codon, oligo_dna)
