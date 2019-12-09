from typing import List, Dict
from collections import OrderedDict
from itertools import chain


def parse_dna(dna_file: str) -> str:
    return [a.rstrip() for a in open(dna_file, "r") if ">" not in a and a != ""][0]


def find_dna_var_poss(var_poss: List[int]) -> List[int]:
    # return list(chain(*([pos * 3 - 2, pos * 3 - 1, pos * 3] for pos in var_poss)))
    all_poss: List[int] = []
    for pos in var_poss:
        all_poss.append(pos * 3 - 2)
        all_poss.append(pos * 3 - 1)
        all_poss.append(pos * 3)
    return all_poss


def expand_dna_var_poss(var_poss: List[int]) -> List[int]:
    return list(chain(*([pos - 2, pos - 1, pos] for pos in var_poss)))


def parse_resfile(in_file: str) -> Dict[int, str]:
    results: Dict[int, str] = OrderedDict()
    for lin in open(in_file, "r"):
        if lin.rstrip() in ["nataa", "start"]:
            continue
        pos = int(lin.split()[0])
        aas = lin.split()[-1]
        results[pos] = aas
    return results
