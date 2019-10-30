from typing import List


def parse_dna(dna_file: str) -> str:
    return [a.rstrip() for a in open(dna_file, "r") if ">" not in a and a != ""][0]


def find_dna_var_poss(var_poss: List[int]) -> List[int]:
    all_poss: List[int] = []
    for pos in var_poss:
        all_poss.append(pos)
        all_poss.append(pos + 1)
        all_poss.append(pos + 2)
    return all_poss