from typing import List

from Bio import SeqFeature


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
