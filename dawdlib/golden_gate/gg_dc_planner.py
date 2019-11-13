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
from collections import OrderedDict

from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gg_path_finder import dijkstra_all_paths, find_compatible_paths
from dawdlib.golden_gate.graph_maker import GraphMaker, SOURCE_NODE, SINK_NODE
from dawdlib.golden_gate.utils import find_dna_var_poss, parse_dna, parse_resfile


from dawdlib.degenerate_dna.codon_selector import CodonSelector

#%%
MIN_OLIGO_LENGTH: int = 20
MAX_OLIGO_LENGTH: int = 80
MIN_CONST_OLIGO_LENGTH: int = 20
MAX_NUM_GATES = 12

gg_data = GGData()
np.random.seed(0)

dna = parse_dna(f"{os.path.dirname(__file__)}/tests/input_921/921_DNA_seq")
resfile = parse_resfile(f"{os.path.dirname(__file__)}/tests/input_921/resfile_921")
aa_var_poss: List[int] = np.array(list(resfile.keys()))
dna_var_poss: List[int] = find_dna_var_poss(aa_var_poss)
const_dna_poss = sorted(set(list(range(1, len(dna)))) - set(dna_var_poss))

#%%
codon_selector = CodonSelector("37762")
aa_pos_n_codons = OrderedDict(
    [
        (pos, codon_selector.optimise_codons(list(aas))["ambiguous_codons"])
        for pos, aas in resfile.items()
    ]
)
dna_pos_n_codons = OrderedDict(
    [(dna_pos, aas) for dna_pos, aas in zip(dna_var_poss[::3], aa_pos_n_codons.values())]
)
#%%
graph_maker = GraphMaker(gg_data=gg_data)
d_graph, src, snk = graph_maker.make_grpah(
    dna=dna,
    const_poss=const_dna_poss,
    var_poss=dna_var_poss,
    dna_pos_n_codons=dna_pos_n_codons,
    min_oligo_length=MIN_OLIGO_LENGTH,
    max_oligo_length=MAX_OLIGO_LENGTH,
    min_const_oligo_length=MIN_CONST_OLIGO_LENGTH,
)
#%%
# sources = [n for n in d_graph.nodes if n.index < dna_var_poss[0]]
# sinks = [n for n in d_graph.nodes if n.index > dna_var_poss[-1]]

#%%
# best_path = next(dijkstra_all_paths(d_graph, sources, sinks, MAX_NUM_GATES, gg_data))
best_paths = nx.all_shortest_paths(d_graph, src, snk, "weight")
best_path: List[Gate] = []
for path in best_paths:
    if not gg_data.gate_set_has_off_target([n.bps for n in path[1:-1]], threshold=1000):
        best_path = path
        break
print(best_path)
#%%
path_edges = [(n1, n2) for n1, n2 in zip(best_path[:-1], best_path[1:])]

fig = plt.figure(figsize=(15, 10))  # , facecolor="w")
ax = plt.axes()
ax.set_facecolor("w")

node_poss = {n: (n.index, np.random.random()) for n in d_graph.nodes}
nx.draw_networkx_nodes(
    d_graph,
    node_poss,
    nodelist=[n for n in d_graph.nodes if n not in best_path],
    node_color="grey",
    node_shape="_",
)
nx.draw_networkx_nodes(
    d_graph,
    node_poss,
    nodelist=best_path,
    node_color="b",
    node_shape=" ",
)
nx.draw_networkx_labels(
    d_graph, node_poss, labels={n: n.bps for n in d_graph.nodes if n in best_path}
)
nx.draw_networkx_edges(
    d_graph, node_poss, edgelist=path_edges, edge_color="b", arrows=True
)

for v_p in dna_var_poss:
    plt.axvline(x=v_p, color="k", linewidth=1)
plt.show()
print(best_path)

# %%


# %%


# %%
