from typing import List

import pandas as pd
from Bio import SeqFeature

from dawdlib.golden_gate.gate import Gate


def create_dc_features(deg_parsed_df):
    features: List[SeqFeature] = []
    for i, row in deg_parsed_df.iterrows():
        num = 0
        for col in [
            c for c in deg_parsed_df.columns if c.startswith("AMBIGUOUS_CODONS")
        ]:
            if str(row[col]) != "nan":
                num += 1
        name = f"{''.join(row['ENCODED_AAS'])}={num}"
        ftr = SeqFeature.SeqFeature(
            SeqFeature.FeatureLocation(row["DNA_POS"] - 1, row["DNA_POS"] + 2, 1),
            type=name,
            id=name,
        )
        features.append(ftr)
    return features


def create_path_features(chosen_path):
    features: List[SeqFeature] = []
    for gate in chosen_path[1:-1]:
        name = f"G.{gate.idx}.{gate.bps}"
        ftr = SeqFeature.SeqFeature(
            SeqFeature.FeatureLocation(gate.idx, gate.idx + 4, 1), type=name, id=name
        )
        features.append(ftr)
    return features


def gate_path_to_df(gate_path):
    gate_df = pd.DataFrame(data=gate_path)
    return gate_df


def df_to_gate_path(gate_df):
    gate_path: List[Gate] = []
    for _, row in gate_df.iterrows():
        gate = Gate(
            idx=row["idx"],
            bps=row["bps"],
            req_primer=row["req_primer"],
            syn_mut=row["syn_mut"],
        )
        gate_path.append(gate)
    return gate_path
