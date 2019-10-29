"""
a notebook for finding the best gates and degenerate codons for a given library
"""
import os

# %%
import sys
from typing import List, Tuple

import networkx as nx

from src.gate_data import GGData

# %% [markdown]
# setup all input

# %%
# W_PATH = ""
# dna_seq = parse_dnaseqs()
# resfile = ResFile
MIN_OLIGO_LENGTH: int = 20
MAX_OLIGO_LENGTH: int = 80
# %%


def parse_dna(dna_file: str) -> str:
    return [a.rstrip() for a in open(dna_file, "r") if ">" not in a and a != ""][0]


def find_dna_des_poss(des_poss: List[int]) -> List[int]:
    all_poss: List[int] = []
    for pos in des_poss:
        all_poss.append(pos)
        all_poss.append(pos + 1)
        all_poss.append(pos + 2)
    return all_poss


# %%
gg_data = GGData()
d_graph = nx.DiGraph()
dna = parse_dna(f"{os.path.dirname(__file__)}/../tests/input_921/921_DNA_seq")
des_poss: List[int] = [30, 31, 33, 56, 57, 59, 99, 101, 102, 104]
dna_des_poss: List[int] = find_dna_des_poss(des_poss)
acceptable_dna_poss = [p for p in range(len(dna)) if p not in dna_des_poss]
# %% [markdown]
# this is where the degenerate codon stuff comes in

# %%


def get_all_self_binding_gates(gg_dsta: GGData) -> List[str]:
    all_gates = gg_data.lig_df.columns
    return [g for g in all_gates if gg_data.score_gate1_rc_gate2(g, g) > 2000]


acceptable_fcw = get_all_self_binding_gates(gg_data)
print(acceptable_fcw, len(acceptable_fcw))
# %%
# find all acceptable gates in the DNA sequence


def is_node(
    fcw: str, ind: int, acceptable_fcw: List[str], acceptable_dna_poss: List[int]
) -> bool:
    if fcw not in acceptable_fcw:
        return False
    if any([p not in acceptable_dna_poss for p in range(ind, ind + 4)]):
        return False
    return True


for ind in range(len(dna) - 3):
    fcw = dna[ind : ind + 4]
    if is_node(fcw, ind, acceptable_fcw, acceptable_dna_poss):
        print(fcw)
        d_graph.add_node((ind, fcw))

# %%
# add all acceptable edges
# TODO add constant region compatability !!!!!!!!!!!


def is_edge(
    nd1: Tuple[int, str], nd2: Tuple[int, str], dna_des_poss: List[int], gg_data: GGData
) -> bool:
    if nd1[0] >= nd2[0]:
        return False

    # segments with variable positions between them are required to be at a certain length
    if any([nd1[0] < p < nd2[0] for p in dna_des_poss]):
        if not (MIN_OLIGO_LENGTH < nd2[0] - nd1[0] < MAX_OLIGO_LENGTH):
            return False
    # if not any([nd1[0] < p < nd2[0] for p in dna_des_poss]):
    #     return False
    if gg_data.gates_all_scores(nd1[1], nd2[1]) > 1000:
        return False
    return True


for nd1 in d_graph.nodes:
    for nd2 in d_graph.nodes:
        if is_edge(nd1, nd2, dna_des_poss, gg_data):
            print(nd1, nd2)
            d_graph.add_edge(nd1, nd2)


# %%
src = (-1, "SRC")
snk = (len(dna) + 1, "SNK")
d_graph.add_node(src)
d_graph.add_node(snk)
for nd1 in d_graph.nodes:
    # if theres a variable position before nd1, conncet it so source only if it is short enough
    # otherwise connect it anyway (constant segment)
    if any([-1 < p < nd1[0] for p in dna_des_poss]):
        if MIN_OLIGO_LENGTH < nd1[0] < MAX_OLIGO_LENGTH:
            d_graph.add_edge(src, nd1)
    else:
        d_graph.add_edge(src, nd1)
    # if theres a variable position between nd1 and sink, conncet it so sink only if it is short enough
    # otherwise connect it anyway (constant segment)
    if any([nd1[0] < p for p in dna_des_poss]):
        if MIN_OLIGO_LENGTH < len(dna) - nd1[0] < MAX_OLIGO_LENGTH:
            d_graph.add_edge(nd1, snk)
    else:
        d_graph.add_edge(nd1, snk)
