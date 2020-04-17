from collections import OrderedDict
from itertools import chain, product
from typing import Dict, Generator, Iterator, List, NamedTuple, Set, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values as adv

from dawdlib.golden_gate.constants import CONST_COST, OLIGO_PREFIX, OLIGO_SUFFIX
from dawdlib.golden_gate.gate import Gate, SynMut


class OligoTableEntry(NamedTuple):
    wt: bool
    const: bool
    full_oligo_dna: str
    gate1: Gate
    gate2: Gate
    gate_gate_dist: int
    name: str
    oligo_codons: List[Tuple[int, str]]
    oligo_dna: str


class Requirements(NamedTuple):
    gg_temp: int
    gg_hours: int
    min_oligo_length: int
    max_oligo_length: int
    min_const_oligo_length: int
    min_efficiency: float
    min_fidelity: float
    oligo_prefix: str  # = OLIGO_PREFIX
    oligo_suffix: str  # = OLIGO_SUFFIX
    const_cost: int  # = CONST_COST
    oligo_addition: int = 0
    re_calc_lengths: bool = True
    filter_gc_overhangs: bool = True


def RequirementsFactory(*args, **kwargs) -> Requirements:
    reqs = Requirements(*args, **kwargs)
    if reqs.re_calc_lengths:
        oligo_addition = len(reqs.oligo_prefix) + len(reqs.oligo_suffix)
        return reqs._replace(
            oligo_addition=oligo_addition,
            max_oligo_length=reqs.max_oligo_length - oligo_addition,
        )


# class Requirements:
#     def __init__(
#         self,
#         gg_temp: int,
#         gg_hours: int,
#         min_oligo_length: int,
#         max_oligo_length: int,
#         min_const_oligo_length: int,
#         min_efficiency: float,
#         min_fidelity: float,
#         oligo_prefix: str = OLIGO_PREFIX,
#         oligo_suffix: str = OLIGO_SUFFIX,
#         re_calc_lengths: bool = True,
#         const_cost=CONST_COST,
#         filter_gc_overhangs: bool = True,
#     ):
#         self.gg_temp = gg_temp
#         self.gg_hours = gg_hours
#         self.min_oligo_length = min_oligo_length
#         self.max_oligo_length = max_oligo_length
#         self.min_const_oligo_length = min_const_oligo_length
#         self.min_efficiency = min_efficiency
#         self.min_fidelity = min_fidelity
#         self.oligo_addition = len(oligo_prefix) + len(oligo_suffix)
#         self.const_cost = const_cost
#         self.filter_gc_overhangs = filter_gc_overhangs
#         if re_calc_lengths:
#             self.max_oligo_length -= self.oligo_addition


def parse_dna(dna_file: str, frmt="fasta") -> str:
    record = SeqIO.read(dna_file, frmt)
    return str(record.seq)


def convert_to_embl(dna_file: str, embl_file: str, frmt="fasta") -> None:
    record = SeqIO.read(dna_file, frmt)
    SeqIO.write(record, embl_file, "embl")


def ambiguous_dna_unambiguous(dna: str) -> Iterator[str]:
    """

    Args:
        dna: An ambiguous dna string

    Returns:
        An iterator over unambiguous DNAs encoded by the input DNA

    """
    dnas = [tuple(adv[char]) for char in dna]
    return map(lambda x: "".join(x), product(*dnas))


def find_dna_var_poss(var_poss: List[int]) -> List[int]:
    all_poss: List[int] = []
    for pos in var_poss:
        all_poss.append(pos * 3 - 2)
        all_poss.append(pos * 3 - 1)
        all_poss.append(pos * 3)
    return all_poss


def expand_dna_var_poss(var_poss: List[int]) -> List[int]:
    return list(chain(*([pos - 1, pos, pos + 1] for pos in var_poss)))


def parse_resfile(in_file: str) -> Dict[int, str]:
    results: Dict[int, str] = OrderedDict()
    for lin in open(in_file):
        if lin.rstrip() in ["nataa", "start"]:
            continue
        pos = int(lin.split()[0])
        aas = lin.split()[-1]
        results[pos] = aas
    return results


def codon_table() -> Dict[str, str]:
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    cdn_table = dict(zip(codons, amino_acids))
    return cdn_table


def amino_acid_codon_table() -> Dict[str, Set[str]]:
    table: Dict[str, Set[str]] = {}
    for key, val in codon_table().items():
        table.setdefault(val, set()).add(key)
    return table


def syn_muts(dna: str) -> Generator[Gate, None, None]:
    # synonymous mutations
    # gate section is from "|" to "|", the codons are seperated by ":"
    # |---:-|-- happens when ind %3 = 0
    # -|--:--|- happens when ind %3 = 1
    # --|-:---| happens when ind %3 = 2
    slices = [slice(0, 4), slice(1, 5), slice(2, 6)]
    cdn_trans = codon_table()
    aa_trans = amino_acid_codon_table()
    for idx in range(0, len(dna) - 3, 3):
        aa1 = cdn_trans[dna[idx : idx + 3]]
        cdn1s = aa_trans[aa1]
        aa2 = cdn_trans[dna[idx + 3 : idx + 6]]
        cdn2s = aa_trans[aa2]

        for bps_slice in slices:
            for cdn1, cdn2 in product(cdn1s, cdn2s):
                yield Gate(
                    idx,
                    (cdn1 + cdn2)[bps_slice],
                    req_primer=True,
                    syn_mut=(SynMut(idx, aa1, cdn1), SynMut(idx + 4, aa2, cdn2)),
                )


def gate_list_df(gate_path: List[Gate]) -> pd.DataFrame:
    return pd.DataFrame.from_records(gate_path, columns=gate_path[0]._fields)


def gate_df_list(gate_df: pd.DataFrame) -> List[Gate]:
    gate_path: List[Gate] = []
    for row in gate_df.itertuples(index=False):
        gate_path.append(Gate(**row._asdict()))
    return gate_path
