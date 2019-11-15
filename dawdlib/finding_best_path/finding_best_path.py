#%%
import pickle
from typing import Dict, Generator, List, Set, Tuple
import random

import networkx as nx
import numpy as np
import pandas as pd

from dawdlib.gg_dc_combine.gg_dc_combine import parse_degenerate_codon_csv
from dawdlib.golden_gate.gate import Gate

# %%
di_graph = nx.read_gpickle(
    "/home/labs/fleishman/jonathaw/gg_n_dc/dawdlib/dawdlib/finding_best_path/tests/921_var_segments_graph.pickle"
)
dc_df = parse_degenerate_codon_csv(
    "/home/labs/fleishman/jonathaw/gg_n_dc/dawdlib/dawdlib/finding_best_path/tests/921.csv"
)
edge_data = pickle.load(
    open(
        "/home/labs/fleishman/jonathaw/gg_n_dc/dawdlib/dawdlib/finding_best_path/tests/921_var_edge_props.pickle",
        "rb",
    )
)

# %%
def find_end_nodes(
    di_graph: nx.DiGraph, dc_df: pd.DataFrame
):  # -> Generator[Tuple[Gate]]:
    last_v_pos: int = dc_df.DNA_POS.max()
    for node in di_graph.nodes:
        if node[1].idx > last_v_pos:
            yield node


# %%
end_nodes = find_end_nodes(di_graph, dc_df)

# %%


def recursively_assign_gates(di_graph: nx.DiGraph, first_v_pos: int, last_v_pos: int):
    count_children = 0
    count_sinks = 0
    nx.set_node_attributes(di_graph, list(), "bpss")
    visited = set()
    stack = [n for n, d in di_graph.in_degree() if d == 0]
    while stack:
        node = stack[-1]

        # reached end of graph
        if di_graph.out_degree[node] == 0:
            count_sinks += 1
            bpss = [node[0].bps, node[1].bps]
            nx.set_node_attributes(di_graph, {node: {"bpss": bpss}})
            stack.pop()
            visited.add(node)
            continue

        # node is visited if all it's children are visited
        if node in visited:
            bpss = [node[0].bps, node[1].bps]
            for child in di_graph.successors(node):
                bpss.extend([child[0].bps])
                bpss.extend([child[1].bps])
            nx.set_node_attributes(di_graph, {node: {"bpss": bpss}})
            stack.pop()
            continue

        # if node is not visited, than the children need to be visited
        for child in di_graph.successors(node):
            if child in visited:
                continue
            count_children += 1
            stack.append(child)
        visited.add(node)
    return count_children, count_sinks


# %%
count_children, count_sinks = recursively_assign_gates(
    di_graph, dc_df.DNA_POS.min(), dc_df.DNA_POS.max()
)
print(count_children)
print(count_sinks)

#%%
def find_v_pos_num_codons(dc_df: pd.DataFrame) -> Dict[int, int]:
    amb_cod_columns = [c for c in dc_df.columns if c.startswith("AMBIGUOUS_CODONS")]
    return {
        p: len(
            [
                a
                for a in list(dc_df.loc[dc_df.DNA_POS == p, amb_cod_columns].values[0])
                if not a is np.nan
            ]
        )
        for p in dc_df.DNA_POS
    }


v_pos_n_codons = find_v_pos_num_codons(dc_df)
print(v_pos_n_codons)
# %%
# def score_path(opt_path: List[Tuple[Gate, Gate]], v_pos_n_codons: Dict[int, int]) -> int:
#     score = 0
#     for node in opt_path:
#         print(node)
#         v_poss = dc_df.loc[dc_df.DNA_POS.between(node[0].idx, node[1].idx)], "DNA_POS"]
#         print(v_poss)
#         score += 0

# print(score_path([], v_pos_n_codons))


#%%
def get_edge_score(rec_path: List[Tuple[Gate, Gate]], node: Tuple[Gate, Gate]) -> int:
    if len(rec_path) == 0:
        return 0
    return edge_data[(rec_path[-1], node)]


def find_all_paths(
    di_graph: nx.DiGraph
) -> Tuple[Dict[int, List[Tuple[Gate, Gate]]], Dict[int, int]]:
    best_score: Dict[int, int] = {l: np.inf for l in range(1, 21)}
    best_path: Dict[int, List[Tuple[Gate, Gate]]] = {l: [] for l in range(1, 21)}

    stack: List[Tuple[Tuple[Gate, Gate], int]] = [
        (n, edge_data[n]["weight"]) for n, d in di_graph.in_degree() if d == 0
    ]
    rec_path: List[Tuple[Gate, Gate]] = []
    visited: Set[Tuple[Gate, Gate]] = set()

    i = 0
    while stack:
        node, upto_node_score = stack[-1]
        if node in visited:
            visited.remove(node)
            rec_path.pop()
            stack.pop()
            continue
        with_node_score = upto_node_score + edge_data[node]["weight"]
        rec_path.append(node)
        if with_node_score > best_score[len(rec_path)]:
            stack.pop()
            rec_path.pop()
            continue

        # reached end of graph
        if di_graph.out_degree[node] == 0:
            stack.pop()
            if best_score[len(rec_path)] > with_node_score:
                best_score[len(rec_path)] = with_node_score
                best_path[len(rec_path)] = rec_path[:]
            rec_path.pop()
            continue

        # if node is not visited, than the children need to be visited
        for child in di_graph.successors(node):
            stack.append((child, with_node_score))
        visited.add(node)
    return best_path, best_score


# %%
n = 50
keep_nodes = random.sample([n for n, d in di_graph.in_degree() if d == 0], n)
keep_nodes += random.sample(
    [n for n, d in di_graph.in_degree() if d != 0 and di_graph.out_degree(n) != 0], n
)
keep_nodes += random.sample([n for n, d in di_graph.out_degree() if d == 0], n)
print(keep_nodes)
sub_graph = di_graph.subgraph(keep_nodes)

print(len(sub_graph.edges))
#%%
best_path, best_score = find_all_paths(di_graph)
print(best_path)
print(best_score)

# %%
