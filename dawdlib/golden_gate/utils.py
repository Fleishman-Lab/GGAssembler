from collections import OrderedDict
from itertools import chain, product
from typing import Dict, Generator, Iterator, List, NamedTuple, Set, Tuple

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


class Requirements:
    def __init__(
        self,
        min_oligo_length: int,
        max_oligo_length: int,
        min_const_oligo_length: int,
        gate_self_binding_min: int = 2000,
        gate_crosstalk_max: int = 1000,
        oligo_prefix: str = OLIGO_PREFIX,
        oligo_suffix: str = OLIGO_SUFFIX,
        re_calc_lengths: bool = True,
        const_cost=CONST_COST,
    ):
        self.min_oligo_length = min_oligo_length
        self.max_oligo_length = max_oligo_length
        self.min_const_oligo_length = min_const_oligo_length
        self.gate_self_binding_min = gate_self_binding_min
        self.gate_crosstalk_max = gate_crosstalk_max
        self.oligo_addition = len(oligo_prefix) + len(oligo_suffix)
        self.const_cost = const_cost
        if re_calc_lengths:
            self.min_oligo_length -= self.oligo_addition
            self.max_oligo_length -= self.oligo_addition


def parse_dna(dna_file: str, frmt="fasta") -> str:
    record = SeqIO.read(dna_file, frmt)
    return str(record.seq)


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
    # return list(chain(*([pos * 3 - 2, pos * 3 - 1, pos * 3] for pos in var_poss)))
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
                    True,
                    (SynMut(idx, aa1, cdn1), SynMut(idx + 4, aa2, cdn2)),
                )
