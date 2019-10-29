"""
a notebook for finding the best gates and degenerate codons for a given library
"""
#%%
import os

import sys
from typing import List, Tuple


from dawdlib.gate_data.gate_data import GGData

# %% [markdown]
# setup all input

# %%
# W_PATH = ""
# dna_seq = parse_dnaseqs()
# resfile = ResFile
MIN_OLIGO_LENGTH: int = 20
MAX_OLIGO_LENGTH: int = 80
MIN_CONST_OLIGO_LENGTH: int = 20
# %%


def parse_dna(dna_file: str) -> str:
    return [a.rstrip() for a in open(dna_file, "r") if ">" not in a and a != ""][0]



# %%
dna = parse_dna(f"{os.path.dirname(__file__)}/../tests/input_921/921_DNA_seq")
aa_var_poss: List[int] = [30, 31, 33, 56, 57, 59, 99, 101, 102, 104]
dna_var_poss: List[int] = find_dna_des_poss(aa_var_poss)
const_dna_poss = [p for p in range(len(dna)) if p not in dna_var_poss]
# %% [markdown]
# this is where the degenerate codon stuff comes in

# %%


acceptable_fcws = gg_data.get_all_self_binding_gates()
print(acceptable_fcws, len(acceptable_fcws))
# %%
# find all acceptable gates in the DNA sequence


def is_node(
    fcw: str, ind: int, acceptable_fcws: List[str], const_dna_poss: List[int]
) -> bool:
    if fcw not in acceptable_fcws:
        return False
    if any([p not in const_dna_poss for p in range(ind, ind + 4)]):
        return False
    return True


for ind in range(len(dna) - 3):
    fcw = dna[ind : ind + 4]
    if is_node(fcw, ind, acceptable_fcws, const_dna_poss):
        print(fcw)
        d_graph.add_node((ind, fcw))

# %%
# add all acceptable edges
# TODO add constant region compatability !!!!!!!!!!!

