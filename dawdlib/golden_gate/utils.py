from typing import List, Dict, Set, Generator
from collections import OrderedDict
from itertools import chain, product
from Bio import SeqIO
from .gate import SynMut, Gate


def parse_dna(dna_file: str, frmt = 'fasta') -> str:
    record = SeqIO.read(dna_file, frmt)
    return str(record.seq)


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
    for lin in open(in_file, "r"):
        if lin.rstrip() in ["nataa", "start"]:
            continue
        pos = int(lin.split()[0])
        aas = lin.split()[-1]
        results[pos] = aas
    return results


def codon_table() -> Dict[str, str]:
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
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
        aa1 = cdn_trans[dna[idx:idx + 3]]
        cdn1s = aa_trans[aa1]
        aa2 = cdn_trans[dna[idx+3:idx + 6]]
        cdn2s = aa_trans[aa2]

        for bps_slice in slices:
            for cdn1, cdn2 in product(cdn1s, cdn2s):
                yield Gate(
                    idx,
                    (cdn1+cdn2)[bps_slice],
                    True,
                    (
                        SynMut(idx, aa1, cdn1),
                        SynMut(idx+4, aa2, cdn2)
                    )
                )
