"""
a notebook for finding the best gates and degenerate codons for a given library
"""
#%%
import os
import sys
from typing import List, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gg_path_finder import dijkstra_all_paths, find_compatible_paths
from dawdlib.golden_gate.graph_maker import GraphMaker
from dawdlib.golden_gate.utils import find_dna_var_poss, parse_dna

#%%
MIN_OLIGO_LENGTH: int = 20
MAX_OLIGO_LENGTH: int = 80
MIN_CONST_OLIGO_LENGTH: int = 20
MAX_NUM_GATES = 12

gg_data = GGData()
np.random.seed(0)

dna = parse_dna(f"{os.path.dirname(__file__)}/tests/input_921/921_DNA_seq")
aa_var_poss: List[int] = [30, 31, 33, 56, 57, 59, 99, 101, 102, 104]
dna_var_poss: List[int] = find_dna_var_poss(aa_var_poss)
const_dna_poss = [p for p in range(len(dna)) if p not in dna_var_poss]

#%%
graph_maker = GraphMaker(gg_data=gg_data)
d_graph: nx.DiGraph = graph_maker.make_grpah(
    dna=dna,
    const_poss=const_dna_poss,
    var_poss=dna_var_poss,
    min_oligo_length=MIN_OLIGO_LENGTH,
    max_oligo_length=MAX_OLIGO_LENGTH,
    min_const_oligo_length=MIN_CONST_OLIGO_LENGTH,
)
#%%
sources = [n for n in d_graph.nodes if n.index < dna_var_poss[0]]
sinks = [n for n in d_graph.nodes if n.index > dna_var_poss[-1]]

#%%
best_path = next(dijkstra_all_paths(d_graph, sources, sinks, MAX_NUM_GATES, gg_data))
print(best_path)
# %% [markdown]
path_edges = [(n1, n2) for n1, n2 in zip(best_path[:-1], best_path[1:])]

fig = plt.figure(figsize=(15, 10))  # , facecolor="w")
ax = plt.axes()
ax.set_facecolor("w")

node_poss = {n: (n.index, np.random.randint(0, 10)) for n in d_graph.nodes}
nx.draw_networkx_nodes(
    d_graph,
    node_poss,
    nodelist=[n for n in d_graph.nodes if n not in best_path],
    node_color="grey",
)
nx.draw_networkx_nodes(
    d_graph,
    node_poss,
    nodelist=best_path,
    node_color="b",
    marker=" ",
)
nx.draw_networkx_labels(d_graph, node_poss, labels={n: n.bps for n in d_graph.nodes if n in best_path})
nx.draw_networkx_edges(
    d_graph, node_poss, edgelist=path_edges, edge_color="b", arrows=True
)

for v_p in dna_var_poss:
    plt.axvline(x=v_p, color="k", linewidth=1)
plt.show()

# %%


# %%
