from collections import OrderedDict
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from dawdlib.degenerate_dna.codon_selector import CodonSelector, PosCodon
from dawdlib.degenerate_dna.structs import TableColNames
from dawdlib.golden_gate.utils import find_dna_var_poss, parse_resfile


def aas_deg_codons(codon_selector: CodonSelector, aas: List[str]) -> PosCodon:
    return codon_selector.optimise_codons(aas)
    # length = len(aas)
    # if length > 15:
    #     return codon_selector.optimise_codons(aas, mode="exact")
    # if length < 4:
    #     return codon_selector.optimise_codons(aas, mode="exact")
    # return codon_selector.optimise_codons(aas, mode="exact")


def resfile_aa_codons(
    codon_selector: CodonSelector, resfile: Dict[int, str]
) -> Dict[int, PosCodon]:
    return OrderedDict(
        [
            (pos, aas_deg_codons(codon_selector, list(aas)))
            for pos, aas in resfile.items()
        ]
    )


def aa_pos_aa_freq(aa_pos_codons: Dict[int, PosCodon]) -> Dict[int, Tuple[List, List]]:
    aa_pos_count = {}
    for pos, codons in aa_pos_codons.items():
        encoded_aas = [
            amino_acid.amino_acid
            for amino_acid in codons.amino_acids
            for _ in amino_acid.codons
        ]
        uniq, counts = np.unique(encoded_aas, return_counts=True)
        aa_pos_count[pos] = (uniq.tolist(), counts.tolist())
    return aa_pos_count


def aa_pos_ambiguous_codons(
    aa_pos_codons: Dict[int, PosCodon]
) -> Dict[int, Tuple[str, ...]]:
    return OrderedDict(
        [(pos, codons.ambiguous_codons) for pos, codons in aa_pos_codons.items()]
    )


def dna_pos_ambiguous_codons(
    aa_pos_deg_codons: Dict[int, Tuple[str, ...]], dna_var_poss: List[int]
) -> Dict[int, Tuple[str, ...]]:
    return OrderedDict(
        (dna_pos, aas)
        for dna_pos, aas in zip(dna_var_poss[::3], aa_pos_deg_codons.values())
    )


def create_deg_table(res_filename: str, codon_selector: CodonSelector) -> pd.DataFrame:
    resfile = parse_resfile(res_filename)
    aa_var_poss: List[int] = list(resfile.keys())
    dna_var_pos: List[int] = find_dna_var_poss(aa_var_poss)

    aa_pos_codons = resfile_aa_codons(codon_selector, resfile)
    aa_pos_deg_codons = aa_pos_ambiguous_codons(aa_pos_codons)
    aa_pos_count = aa_pos_aa_freq(aa_pos_codons)
    dna_pos_n_codons = dna_pos_ambiguous_codons(aa_pos_deg_codons, dna_var_pos)

    aa_pos = aa_pos_deg_codons.keys()
    dna_pos = dna_pos_n_codons.keys()
    ambiguous_codons = aa_pos_deg_codons.values()

    df = pd.DataFrame((aa_pos, dna_pos)).T
    df.columns = [TableColNames.AA_POS.value, TableColNames.DNA_POS.value]
    ambiguous_codons_df = pd.DataFrame(ambiguous_codons)
    ambiguous_codons_df.columns = [
        f"{TableColNames.AMBIGUOUS_CODONS.value}{i+1}"
        for i in ambiguous_codons_df.columns
    ]
    aa_pos_count_df = pd.DataFrame(
        aa_pos_count.values(),
        columns=[TableColNames.ENCODED_AAS.value, TableColNames.ENCODED_COUNT.value],
    )
    df = pd.concat((df, aa_pos_count_df, ambiguous_codons_df), axis=1)
    return df


def generate_deg_csv(
    res_filename: str,
    csv_filename: str,
    organism_id: str = "37762",
    mode: str = "exact",
):
    """
    Default organism is Escherichia coli, ID = "37762"
    Saccharomyces cerevisiae ID = "4932"

    For to search for a different organism
        see: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?
        https://doi.org/10.1093/nar/gkr1178
    """
    cs = CodonSelector(organism_id, mode)
    df: pd.DataFrame = create_deg_table(res_filename, cs)
    df.to_csv(csv_filename, index=False)
