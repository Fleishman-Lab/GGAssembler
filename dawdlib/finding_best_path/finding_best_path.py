#%%
# %%
import cProfile
import pickle
# %%
import pstats
import random
import time
import timeit
from typing import Callable, Dict, Generator, List, Set, Tuple

import networkx as nx
import numpy as np
import pandas as pd

from dawdlib.gg_dc_combine.gg_dc_combine import parse_degenerate_codon_csv
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gate_restriction import create_path_validator

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

#%%
di_graph_4 = nx.read_gpickle(
    "/home/labs/fleishman/jonathaw/gg_n_dc/dawdlib/dawdlib/finding_best_path/tests/921_4_var_segments_graph.pickle"
)
edge_data_4 = pickle.load(
    open(
        "/home/labs/fleishman/jonathaw/gg_n_dc/dawdlib/dawdlib/finding_best_path/tests/921_4_var_edge_props.pickle",
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
def recursively_assign_subgraph_scores(di_graph: nx.DiGraph, edge_data) -> None:
    nx.set_node_attributes(di_graph, 0, "subgraph_score")
    visited = set()
    stack = [n for n, d in di_graph.in_degree() if d == 0]
    while stack:
        node = stack[-1]

        # reached end of graph
        if di_graph.out_degree[node] == 0:
            nx.set_node_attributes(
                di_graph, {node: {"subgraph_score": edge_data[node]["weight"]}}
            )
            stack.pop()
            visited.add(node)
            continue

        # node is visited if all it's children are visited
        if node in visited:
            sg_score = edge_data[node]["weight"]
            for child in di_graph.successors(node):
                sg_score += di_graph.nodes[child]["subgraph_score"]
            nx.set_node_attributes(di_graph, {node: {"subgraph_score": sg_score}})
            stack.pop()
            continue

        # if node is not visited, than the children need to be visited
        for child in di_graph.successors(node):
            if child in visited:
                continue
            stack.append(child)
        visited.add(node)


recursively_assign_subgraph_scores(di_graph_4, edge_data=edge_data_4)
recursively_assign_subgraph_scores(di_graph, edge_data=edge_data)

#%%


def validate_node_path(
    path_validator: Callable, node_path: List[Tuple[Gate, Gate]]
) -> bool:
    gates = []
    for node in node_path:
        gates.append(node[0])
        gates.append(node[1])
    return path_validator(gates)


def find_all_paths(
    di_graph: nx.DiGraph, edge_data: Dict[Tuple[Gate, Gate], int], max_nodes: int = 8
) -> Tuple[Dict[int, List[Tuple[Gate, Gate]]], Dict[int, int]]:
    gg_data = GGData()
    LEARN = 100

    p_valid = create_path_validator(gg_data, 2000, 1000)

    best_score: Dict[int, List[int]] = {l: [np.inf] for l in range(1, max_nodes + 1)}
    best_path: Dict[int, List[Tuple[Gate, Gate]]] = {
        l: [] for l in range(1, max_nodes + 1)
    }

    pre_stack: List[Tuple[int, Tuple[Gate, Gate]]] = sorted(
        [(edge_data[n]["weight"], n) for n, d in di_graph.in_degree() if d == 0]
    )
    stack: List[Tuple[Tuple[Gate, Gate], int]] = [(n[1], 0) for n in pre_stack]
    rec_path: List[Tuple[Gate, Gate]] = []
    visited: Set[Tuple[Gate, Gate]] = set()

    finished: Dict[int, bool] = {}

    while stack:
        node, upto_node_score = stack[-1]
        if node in visited:
            visited.remove(node)
            rec_path.pop()
            stack.pop()
            continue
        with_node_score = upto_node_score + edge_data[node]["weight"]
        rec_path.append(node)
        # reached end of graph
        if di_graph.out_degree[node] == 0:
            stack.pop()
            # print("reached end")
            if best_score[len(rec_path)][-1] > with_node_score:
                best_score[len(rec_path)].append(with_node_score)
                best_path[len(rec_path)] = rec_path[:]
                if len(best_score[len(rec_path)]) > LEARN:
                    finished[len(rec_path)] = all(
                        [
                            sc1 > 0.995 * sc2
                            for sc1, sc2 in zip(
                                best_score[len(rec_path)][-LEARN-1:],
                                best_score[len(rec_path)][-LEARN:-1],
                            )
                        ]
                    )
                    print(best_score[len(rec_path)][-10:], all(list(finished.values())))
                    if all(list(finished.values())):
                        break

            rec_path.pop()
            continue

        # dead-end elimination - if the path up to node has a worse score,
        # or the gates are not comlementary, remove this path
        if (
            len(rec_path) + 1 > max_nodes
            or all(
                [
                    with_node_score > best_score[size_][-1]
                    for size_ in range(len(rec_path) + 1, max_nodes + 1)
                ]
            )
            or not validate_node_path(p_valid, rec_path)
        ):
            stack.pop()
            rec_path.pop()
            continue

        # if node is not visited, than the children need to be visited
        srtd_children = sorted(
            [
                (di_graph.nodes[n]["subgraph_score"], n)
                for n in di_graph.successors(node)
            ]
        )

        for _, child in srtd_children:
            stack.append((child, with_node_score))
        visited.add(node)
    return best_path, {k: v[-1] for k, v in best_score.items()}


# %%
n = 10
keep_nodes = random.sample([n for n, d in di_graph.in_degree() if d == 0], n)
keep_nodes += random.sample(
    [n for n, d in di_graph.in_degree() if d != 0 and di_graph.out_degree(n) != 0], n
)
keep_nodes += random.sample([n for n, d in di_graph.out_degree() if d == 0], n)
# print(keep_nodes)
sub_graph = di_graph.subgraph(keep_nodes)

print(len(sub_graph.edges))
start = timeit.default_timer()
best_paths, best_score = find_all_paths(sub_graph, edge_data)
stop = timeit.default_timer()
# %timeit -n 1 -r 1 find_all_paths(sub_graph, edge_data)
for length in range(1, 13):
    if length in best_paths.keys():
        if best_paths[length]:
            print(length)
            print(best_paths[length])
            print(best_score[length])
print("stop - start", stop - start)
#%%
start = timeit.default_timer()
best_paths, best_score = find_all_paths(di_graph, edge_data)
stop = timeit.default_timer()
for length in range(1, 13):
    if length in best_paths.keys():
        print(length)
        print(best_paths[length])
        print(best_score[length])
print("stop - start", stop - start)
# %%


# %timeit -n 1 -r 1 best_paths, best_score = find_all_paths(di_graph_4, edge_data_4)
start = timeit.default_timer()
best_paths, best_score = find_all_paths(di_graph_4, edge_data_4, max_nodes=13)
stop = timeit.default_timer()
for length in range(1, 13):
    if length in best_paths.keys():
        if best_paths[length]:
            print(length)
            print(best_paths[length])
            print(best_score[length])
print("stop - start", stop - start)

# %%

# %%



def t(a):
    time.sleep(5000)
    print(111111111111111111)
    return False


if True or t("A"):
    print(222222222)


cProfile.run(
    "find_all_paths(di_graph_4, edge_data_4)", "/home/labs/fleishman/jonathaw/aaaa4"
)


p = pstats.Stats("/home/labs/fleishman/jonathaw/aaaa")
p.sort_stats("cumulative").print_stats(100)

# %%
