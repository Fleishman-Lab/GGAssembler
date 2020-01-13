# import typing as tp
# from functools import reduce
# from heapq import heapify, heappop, heappush
# from itertools import combinations, product, count
# from collections import defaultdict
# from enum import Enum, auto
# from functools import lru_cache
#
# import tqdm
#
# import pickle
# import networkx as nx
# import numpy as np
# from dawdlib.golden_gate.gate import Gate
# from dawdlib.golden_gate.gate_data import GGData
# from dawdlib.golden_gate.graph_maker import GraphMaker, make_default_graph
# from networkx.algorithms.approximation import treewidth
#
# import os
# from itertools import chain
# from typing import Callable, Generator, Iterable, List, Dict
# import logging
# import networkx as nx
# import pandas as pd
# from dawdlib.golden_gate.dijkstra import all_shortest_paths
# from dawdlib.golden_gate.gate import Gate, GateSet
# from dawdlib.golden_gate.gate_restriction import create_path_validator
# from dawdlib.golden_gate.gate_data import GGData
# from dawdlib.golden_gate.graph_maker import make_default_graph, GraphMaker
# from dawdlib.degenerate_dna.deg_table import TableColNames
# from dawdlib.golden_gate.utils import parse_dna, expand_dna_var_poss
# from dawdlib.golden_gate.find_gg import deg_table_to_dict
#
# from collections import defaultdict
# from functools import partial
# from itertools import repeat
# from heapq import heappush, heappop, heapify
#
#
# class Bag(tp.NamedTuple):
#     vertices: frozenset
#     id: int
#
#
# class BagType(Enum):
#     Leaf = auto()
#     Introduce = auto()
#     Forget = auto()
#     Join = auto()
#
#
# # class Leaf(Bag):
# #     type = "leaf"
#
#
# # class Introduce(Bag):
# #     type = "introduce"
#
#
# # class Forget(Bag):
# #     type = "forget"
#
#
# # class Join(Bag):
# #     type = "join"
#
#
# def nested_defaultdict(default_factory, depth=1):
#     result = partial(defaultdict, default_factory)
#     for _ in repeat(None, depth - 1):
#         result = partial(defaultdict, result)
#     return result()
#
#
# class GateMinDegreeHeuristic:
#     """ Implements the Minimum Degree heuristic.
#
#     The heuristic chooses the nodes according to their degree
#     (number of neighbours), i.e., first the node with the lowest degree is
#     chosen, then the graph is updated and the corresponding node is
#     removed. Next, a new node with the lowest degree is chosen, and so on.
#     """
#
#     def __init__(self, graph):
#         self._graph = graph
#
#         # nodes that have to be updated in the heap before each iteration
#         self._update_nodes = []
#
#         self._degreeq = []  # a heapq with 3-tuples (degree, gate index, gate)
#
#         # build heap with initial degrees
#         for n in graph:
#             self._degreeq.append((len(graph[n]), n.index, n))
#         heapify(self._degreeq)
#
#     def best_node(self, graph):
#         # update nodes in self._update_nodes
#         for n in self._update_nodes:
#             # insert changed degrees into degreeq
#             heappush(self._degreeq, (len(graph[n]), n.index, n))
#
#         # get the next valid (minimum degree) node
#         while self._degreeq:
#             (min_degree, node_index, elim_node) = heappop(self._degreeq)
#             if elim_node not in graph or len(graph[elim_node]) != min_degree:
#                 # outdated entry in degreeq
#                 continue
#             elif min_degree == len(graph) - 1:
#                 # fully connected: abort condition
#                 return None
#
#             # remember to update nodes in the heap before getting the next node
#             self._update_nodes = graph[elim_node]
#             return elim_node
#
#         # the heap is empty: abort
#         return None
#
#
# def find_suitable_root(tree: nx.Graph, dummy_node):
#     """
#     Not sure if that is the right thing to do
#     But we'll return the bag with the least degree
#     """
#     # return next(iter(nx.topological_sort(tree)))
#     max_degree = np.NINF
#     max_node = None
#     for node, deg in tree.degree():
#         if node.vertices == frozenset(dummy_node):
#             return node
#         if deg > max_degree:
#             max_degree = deg
#             max_node = node
#     return max_node
#
#
# def contract(tree: nx.Graph, u, v) -> None:
#     delta = set(tree.neighbors(v)) - set([u])
#     tree.remove_node(v)
#     for w in delta:
#         tree.add_edge(u, w)
#
#
# def baggify_tree(tree: nx.Graph) -> tp.Tuple[nx.Graph, tp.Dict, int]:
#     bag_tree = nx.Graph()
#     nodes_dict = {}
#     for i, node in enumerate(tree.nodes):
#         nodes_dict[node] = Bag(node, i)
#     bag_tree.add_nodes_from(nodes_dict.values())
#     edges = [(nodes_dict[u], nodes_dict[v]) for u, v in tree.edges]
#     bag_tree.add_edges_from(edges)
#     return bag_tree, nodes_dict, i
#
#
# def make_nice(suitable_root, tree: nx.Graph, next_id: int):
#     root = Bag(frozenset(), next_id)
#     next_id += 1
#     tree.add_edge(root, suitable_root)
#
#     stack: tp.List[Bag] = []
#     visited: tp.Set[Bag] = set()
#
#     stack.append(root)
#
#     while stack:
#         # nx.drawing.nx_pydot.write_dot(
#         #     tree, "/home/labs/fleishman/shayho/Downloads/tree.dot"
#         # )
#         v = stack.pop()
#         visited.add(v)
#
#         # get non-visited neighbors
#         neighbors = list(filter(lambda x: x not in visited, tree.neighbors(v)))
#         # compute number of non-visited neighbors
#         k = len(neighbors)
#
#         if k == 0 and len(v.vertices) != 0:
#             # we have no children, but we are not empty too
#             leaf = Bag(frozenset(), next_id)
#             next_id += 1
#             tree.add_edge(v, leaf)
#             stack.append(v)
#         elif k == 1:
#             w = neighbors[-1]
#
#             if v.vertices == w.vertices:
#                 contract(tree, v, w)
#                 stack.append(v)
#                 continue
#
#             # otherwise we have to reduce the distance
#             new_bag_vertices = set(v.vertices)
#             delta = v.vertices.difference(w.vertices)
#             if len(delta) > 0:
#                 new_bag_vertices.remove(next(iter(delta)))
#             else:
#                 delta = w.vertices.difference(v.vertices)
#                 new_bag_vertices.add(next(iter(delta)))
#             new_bag = Bag(frozenset(new_bag_vertices), next_id)
#             next_id += 1
#             tree.remove_edge(v, w)
#             tree.add_edge(v, new_bag)
#             tree.add_edge(new_bag, w)
#             stack.append(new_bag)
#         elif k > 1:
#             left = Bag(frozenset(v.vertices), next_id)
#             next_id += 1
#             right = Bag(frozenset(v.vertices), next_id)
#             next_id += 1
#             neighbors = list(filter(lambda x: x not in visited, tree.neighbors(v)))
#             for w in neighbors:
#                 tree.remove_edge(v, w)
#             tree.add_edge(v, left)
#             tree.add_edge(v, right)
#             for i, w in enumerate(neighbors):
#                 if i < len(neighbors) / 2:
#                     tree.add_edge(left, w)
#                 else:
#                     tree.add_edge(right, w)
#             stack.append(left)
#             stack.append(right)
#
#     return root
#
#
# def classify_bags(tree: nx.Graph, root) -> tp.Dict[Bag, BagType]:
#
#     bag_type: tp.Dict[Bag, str] = {}
#
#     stack: tp.List[Bag] = []
#     visited: tp.Set[Bag] = set()
#     visited.add(root)
#     stack.append(root)
#
#     while stack:
#         u = stack.pop()
#
#         if len(tree[u]) == 2:
#             bag_type[u] = BagType.Join
#         for v in tree[u]:
#             if v in visited:
#                 continue
#             if u not in bag_type.keys():
#                 tmp = u.vertices - v.vertices
#                 if len(tmp) == 1:
#                     bag_type[u] = BagType.Introduce
#                 else:
#                     bag_type[u] = BagType.Forget
#             visited.add(v)
#             stack.append(v)
#         if u not in bag_type.keys():
#             bag_type[u] = BagType.Leaf
#
#     return bag_type
#
#
# def powerset(iterable: tp.Iterable) -> Generator[frozenset, None, None]:
#     "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
#     s = frozenset(iterable)
#     return map(
#         frozenset, chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))
#     )
#
#
# def is_independent_set(G: nx.Graph, nodes: tp.Iterable) -> bool:
#     return len(G.subgraph(nodes).edges) == 0
#
#
# def segments_to_bps(segments: tp.List[tp.Tuple[Gate, Gate]]) -> tp.FrozenSet[str]:
#     acc = set()
#     for segment in segments:
#         for gate in segment:
#             acc.add(gate.bps)
#     return frozenset(acc)
#
#
# def AlgTDCIG(
#     RGraph: nx.Graph,
#     Rtree: nx.Graph,
#     root,
#     IG: nx.Graph,
#     bps_segs: Dict[str, tp.Tuple[Gate, Gate]],
#     bag_type: tp.Dict[Bag, str],
#     p: tp.Callable[[tp.Tuple[Gate, Gate]], int],
#     w: tp.Callable[[tp.Tuple[Gate, Gate]], int],
#     min_P: int,
#     max_P: int,
#     C: int,
# ):
#     """[summary]
#
#     Args:
#         G (nx.Graph): DNA Graph - vertices are Gates edges are DNA segments with cost
#         RGraph (nx.Graph): Golden gate restriction graph - vertices are 4 bps edges are incompatibility
#         Rtree (nx.Graph): Golden gate restriction graph tree decomposition
#         root ([type]): restriction graph root
#         IG (nx.Graph): DNA segment graph - vertices are variable segments, edges connect vertices if the segments overlap
#         bps_gates (Dict[str, Gate]): 4 bps to Gates tuple (segment) mapping
#         bag_type (tp.Dict[Bag, str]): Mapping from Golden gate restriction graph nodes to their type
#         p (tp.Callable[tp.Tuple[Gate, Gate], int]): Function which calculates the profit (or cost) of each variable segment
#         min_P: (int): minimum profit to consider
#         max_P: (int): maximum profict to consider,
#         C: (int): maximum weight,
#
#     Returns:
#         [type]: [description]
#     """
#
#     def copyf(J: Bag, R: Bag):
#         for S in f[J].keys():
#             for T in f[J][S].keys():
#                 for d in f[J][S][T].keys():
#                     f[R][S][T][d] = f[J][S][T][d]
#
#     def leaf(R: Bag):
#         v = R.vertices
#         if len(v) != 1:
#             raise ValueError("Vertices count is larger than 1", R)
#         for u_seg in bps_segs[next(iter(v))]:
#             u = frozenset([u_seg])
#             f[R][v][u][p(u_seg)] = (w(u_seg), u)
#         f[R][frozenset()][frozenset()][0] = (0, frozenset())
#
#     def introduce(R: Bag):
#         children: tp.List[Bag] = Rtree[R]
#         if len(children) != 1:
#             raise ValueError("Children count is larger than 1 parent: %s", R)
#         J = next(iter(children))
#         v = R.vertices - J.vertices
#         v_bps = next(iter(v))
#
#         copyf(J, R)
#
#         S_keys = list(powerset(J.vertices))
#         for S in S_keys:
#             if not is_independent_set(RGraph, S | v):
#                 continue
#             for T in f[J][S].keys():
#                 for u_seg in bps_segs[v_bps]:
#                     u = frozenset([u_seg])
#                     if (
#                         len(segments_to_bps(u) & segments_to_bps(T)) == 0
#                         and is_independent_set(IG, T | u)
#                         and is_independent_set(RGraph, segments_to_bps(T | u))
#                         and sum(w(x) for x in T | u) <= C
#                     ):
#                         d_keys = list(map(lambda x: x + p(u_seg), f[J][S][T].keys()))
#                         # we want to update only d for which d-p(u) is in d_keys
#                         # this means only d that is in d_keys+p(u)
#                         for d in d_keys:
#                             f[R][S | v][T | u][d] = (
#                                 w(u_seg) + f[J][S][T][d - p(u_seg)][0],
#                                 u | f[J][S][T][d - p(u_seg)][1],
#                             )
#
#         # for d in range(0, max_P+1):
#         #     S_keys = list(powerset(J.vertices))
#         #     for S in S_keys:
#         #         if not is_independent_set(RGraph, S | v):
#         #             continue
#         #         # H_S = [bps_segs[s] for s in S]
#         #         # # TODO: FIX THIS - there is no function independent_sets
#         #         # # must implement a function which iterates
#         #         # # over independent sets of IG.subgraph(H_S)
#         #         # for T in IG.subgraph(H_S).independent_sets():
#         #         T_keys = list(f[J][d][S].keys())
#         #         for T in T_keys:
#         #             for u_seg in bps_segs[v_bps]:
#         #                 u = frozenset([u_seg])
#         #                 if (
#         #                     T != u
#         #                     and is_independent_set(IG, T | u)
#         #                     and is_independent_set(RGraph, segments_to_bps(T | u))
#         #                     and sum(w(x) for x in T | u) <= C
#         #                 ):
#         #                     f[R][d][S | v][T | u] = (w(u_seg) + f[J][d - p(u_seg)][S][T][0], u | f[J][d - p(u_seg)][S][T][1])
#
#     def forget(R: Bag):
#         children: tp.List[Bag] = Rtree[R]
#         if len(children) != 1:
#             raise ValueError("Children count is larger than 1 parent: %s", R)
#         J = next(iter(children))
#         v = J.vertices - R.vertices
#         v_bps = next(iter(v))
#
#         copyf(J, R)
#
#         S_keys = list(powerset(J.vertices))
#         for S in S_keys:
#             if not is_independent_set(RGraph, S | v):
#                 continue
#             for T in f[J][S].keys():
#                 for u_seg in bps_segs[v_bps]:
#                     u = frozenset([u_seg])
#                     if (
#                         len(segments_to_bps(u) & segments_to_bps(T)) == 0
#                         and is_independent_set(IG, T | u)
#                         and is_independent_set(RGraph, segments_to_bps(T | u))
#                         and sum(w(x) for x in T | u) <= C
#                     ):
#                         d_keys = list(f[J][S][T].keys())
#                         ud_keys = list(f[J][S | v][T | u].keys())
#                         # we want to update only d for which d-p(u) is in d_keys
#                         # this means only d that is in d_keys+p(u)
#                         for ud in ud_keys:
#                             if ud in d_keys:
#                                 if f[J][S | v][T | u][ud] < f[J][S][T][ud]:
#                                     f[R][S][T | u][ud] = f[R][S | v][T | u][ud]
#                             else:
#                                 f[R][S][T | u][ud] = f[R][S | v][T | u][ud]
#
#         # for d in range(0, max_P+1):
#         #     S_keys = list(powerset(J.vertices))
#         #     for S in S_keys:
#         #         if not is_independent_set(RGraph, S | v):
#         #             continue
#         #         # H_S = [bps_segs[s] for s in S]
#         #         # # TODO: FIX THIS - there is no function independent_sets
#         #         # # must implement a function which iterates
#         #         # # over independent sets of IG.subgraph(H_S)
#         #         # for T in IG.subgraph(H_S).independent_sets():
#         #         T_keys = list(f[J][d][S].keys())
#         #         for T in T_keys:
#         #             for u_seg in bps_segs[v_bps]:
#         #                 u = frozenset([u_seg])
#         #                 if (
#         #                     T != u
#         #                     and is_independent_set(IG, T | u)
#         #                     and is_independent_set(RGraph, segments_to_bps(T | u))
#         #                     and sum(w(x) for x in T | u) <= C
#         #                     and f[J][d][S | v][T | u] < f[J][d][S][T]
#         #                 ):
#         #                     f[R][d][S][T | u] = (f[J][d][S | v][T | u][0], u | f[J][d][S | v][T | u][1])
#         #                     # f[R][S][T | u][d] = min([f[J][S][T][d], f[J][S | v][T | u][d]])
#
#     def join(R: Bag):
#         children: tp.List[Bag] = Rtree[R]
#         if len(children) != 2:
#             raise ValueError("Children count is larger than 1 parent: %s", R)
#         children_iter = iter(children)
#         J = next(children_iter)
#         copyf(J, R)
#         J = next(children_iter)
#
#         S_keys = list(powerset(J.vertices))
#         for S in S_keys:
#             TJ_keys = list(f[J][S].keys())
#             TR_keys = list(f[R][S].keys())
#             for TJ, TR in product(TJ_keys, TR_keys):
#                 if (
#                     TJ == TR
#                     or not is_independent_set(IG, TJ | TR)
#                     or not is_independent_set(RGraph, segments_to_bps(TJ | TR))
#                 ):
#                     continue
#                 Jd_keys = list(f[J][S][TJ].keys())
#                 Rd_keys = list(f[R][S][TR].keys())
#                 # We want to update d values for which d-k (k in Jd) is in Rd
#                 # this means d in sum(Rd x Jd) (the cartesian product)
#                 # or there is r in Rd and j in Jd such that d = r + j
#                 comb_keys = defaultdict(set)
#                 for r, j in product(Rd_keys, Jd_keys):
#                     comb_keys[r + j].add((r, j))
#                 for d, rj_set in comb_keys.items():
#                     f[R][S][TR | TJ][d] = min(
#                         [
#                             (
#                                 f[R][S][TR][r][0] + f[J][S][TJ][j][0],
#                                 f[R][S][TR][r][1] | f[J][S][TJ][j][1],
#                             )
#                             for r, j in rj_set
#                         ],
#                         key=lambda x: x[0],
#                     )
#
#         # for d in range(0, max_P+1):
#         #     S_keys = list(powerset(J.vertices))
#         #     for S in S_keys:
#         #         T_keys = list(f[J][d][S].keys())
#         #         for T in T_keys:
#         #             f[R][d][S][T] = min(
#         #                 [
#         #                     (f[R][d-k][S][T][0] + f[J][k][S][T][0], f[R][d-k][S][T][1] | f[J][k][S][T][1])
#         #                     for k in range(0, d+1)
#         #                 ],
#         #                 key=lambda x: x[0]
#         #             )
#
#     stack: tp.List[Bag] = [root]
#     # f is a 4 dimensional matrix:
#     #   1st dimension is Bags (restriction tree vertices)
#     #   2nd dimenstion is possible sets of bps
#     #   3rd dimension is possible sets of variable segments (intervals)
#     #   4th dimension is possible costs
#     # f: tp.Dict[
#     #     Bag,
#     #     tp.Dict[frozenset, tp.Dict[frozenset, tp.Dict[int, tp.Tuple[int, frozenset]]]],
#     # ] = defaultdict(
#     #     lambda: defaultdict(
#     #         lambda: defaultdict(lambda: defaultdict(lambda: (C + 1, frozenset())))
#     #     )
#     # )
#     f: tp.Dict[
#         Bag,
#         tp.Dict[frozenset, tp.Dict[frozenset, tp.Dict[int, tp.Tuple[int, frozenset]]]],
#     ] = nested_defaultdict(dict, depth=3)
#
#     visited: tp.Set[Bag] = set()
#
#     while stack:
#         R = stack[-1]
#         if bag_type[R] == BagType.Leaf:
#             visited.add(R)
#             stack.pop()
#             leaf(R)
#         else:
#             if R in visited:
#                 # Finished processing R's children its time to process R
#                 # If R is Introduce
#                 if bag_type[R] == BagType.Introduce:
#                     introduce(R)
#                 # If R is forget
#                 elif bag_type[R] == BagType.Forget:
#                     forget(R)
#                     pass
#                 # if R is Join
#                 elif bag_type[R] == BagType.Join:
#                     join(R)
#                     pass
#                 stack.pop()
#                 pass
#             else:
#                 # R is not leaf - push R's children to the stack and mark it
#                 # as visited
#                 stack.extend(list(Rtree[R]))
#                 visited.add(R)
#     return f
#
#
# def AlgTDCIG2(
#     RGraph: nx.Graph,
#     Rtree: nx.Graph,
#     root,
#     IG: nx.Graph,
#     bps_segs: tp.Dict[str, tp.Set[tp.Tuple[Gate, Gate]]],
#     bag_type: tp.Dict[Bag, str],
#     p: tp.Callable[[tp.Tuple[Gate, Gate]], int],
#     w: tp.Callable[[tp.Tuple[Gate, Gate]], int],
#     min_P: int,
#     max_P: int,
#     C: int,
# ):
#     """[summary]
#
#     Args:
#         G (nx.Graph): DNA Graph - vertices are Gates edges are DNA segments with cost
#         RGraph (nx.Graph): Golden gate restriction graph - vertices are 4 bps edges are incompatibility
#         Rtree (nx.Graph): Golden gate restriction graph tree decomposition
#         root ([type]): restriction graph root
#         IG (nx.Graph): DNA segment graph - vertices are variable segments, edges connect vertices if the segments overlap
#         bps_gates (Dict[str, Gate]): 4 bps to Gates tuple (segment) mapping
#         bag_type (tp.Dict[Bag, str]): Mapping from Golden gate restriction graph nodes to their type
#         p (tp.Callable[tp.Tuple[Gate, Gate], int]): Function which calculates the profit (or cost) of each variable segment
#         min_P: (int): minimum profit to consider
#         max_P: (int): maximum profict to consider,
#         C: (int): maximum weight,
#
#     Returns:
#         [type]: [description]
#     """
#
#     def copyf(J: Bag, R: Bag):
#         for S in f[J].keys():
#             copyfS(J, R, S)
#
#     def copyfS(J: Bag, R: Bag, S):
#         for T in f[J][S].keys():
#             for d in f[J][S][T].keys():
#                 f[R][S][T][d] = f[J][S][T][d]
#
#     def leaf(R: Bag):
#         v = R.vertices
#         if len(v) != 1:
#             raise ValueError("Vertices count is larger than 1 ", R)
#         v_bps = next(iter(v))
#         while bps_segs[v_bps]:
#             u_seg = bps_segs[v_bps].pop()
#             u = frozenset([u_seg])
#             other_bps = segments_to_bps(u) - v
#             bps_segs[other_bps].discard(u_seg)
#             f[R][v][u][p(u_seg)] = (w(u_seg), u)
#         f[R][frozenset()][frozenset()][0] = (0, frozenset())
#
#     def introduce(R: Bag):
#         children: tp.List[Bag] = Rtree[R]
#         if len(children) != 1:
#             raise ValueError("Children count is larger than 1 ", R)
#         J = next(iter(children))
#         v = R.vertices - J.vertices
#         v_bps = next(iter(v))
#
#         copyf(J, R)
#
#         S_keys = list(powerset(J.vertices))
#         for S in S_keys:
#             if not is_independent_set(RGraph, S | v):
#                 continue
#             for T in f[J][S].keys():
#                 while bps_segs[v_bps]:
#                     u_seg = bps_segs[v_bps].pop()
#                     u = frozenset([u_seg])
#                     other_bps = segments_to_bps(u) - v
#                     bps_segs[other_bps].discard(u_seg)
#                     if (
#                         len(segments_to_bps(u) & segments_to_bps(T)) == 0
#                         and is_independent_set(IG, T | u)
#                         and is_independent_set(RGraph, segments_to_bps(T | u))
#                         and sum(w(x) for x in T | u) <= C
#                     ):
#                         d_keys = list(map(lambda x: x + p(u_seg), f[J][S][T].keys()))
#                         # we want to update only d for which d-p(u) is in d_keys
#                         # this means only d that is in d_keys+p(u)
#                         for d in d_keys:
#                             f[R][S | v][T | u][d] = (
#                                 w(u_seg) + f[J][S][T][d - p(u_seg)][0],
#                                 u | f[J][S][T][d - p(u_seg)][1],
#                             )
#
#     def forget(R: Bag):
#         children: tp.List[Bag] = Rtree[R]
#         if len(children) != 1:
#             raise ValueError("Children count is larger than 1 ", R)
#         J = next(iter(children))
#         v = J.vertices - R.vertices
#
#         # copyf(J, R)
#
#         S_keys = list(powerset(J.vertices))
#         for S in S_keys:
#             if S == v or not is_independent_set(RGraph, S | v):
#                 copyfS(J, R, S)
#                 continue
#             for TS, Tv in product(f[J][S].keys(), f[J][S | v].keys()):
#                 d_keys = set(f[J][S][TS].keys())
#                 ud_keys = set(f[J][S | v][Tv].keys())
#                 for d in d_keys & ud_keys:
#                     if f[J][S | v][Tv][d] < f[J][S][TS][d]:
#                         f[R][S][Tv][d] = f[J][S | v][Tv][d]
#                 for d in ud_keys - d_keys:
#                     f[R][S][Tv][d] = f[J][S | v][Tv][d]
#
#     def join(R: Bag):
#         children: tp.List[Bag] = Rtree[R]
#         if len(children) != 2:
#             raise ValueError("Children count is larger than 1 ", R)
#         children_iter = iter(children)
#         J = next(children_iter)
#         copyf(J, R)
#         J = next(children_iter)
#
#         S_keys = list(powerset(J.vertices))
#         for S in S_keys:
#             TJ_keys = list(f[J][S].keys())
#             TR_keys = list(f[R][S].keys())
#             for TJ, TR in product(TJ_keys, TR_keys):
#                 if (
#                     TJ == TR
#                     or not is_independent_set(IG, TJ | TR)
#                     or not is_independent_set(RGraph, segments_to_bps(TJ | TR))
#                 ):
#                     continue
#                 Jd_keys = list(f[J][S][TJ].keys())
#                 Rd_keys = list(f[R][S][TR].keys())
#                 # We want to update d values for which d-k (k in Jd) is in Rd
#                 # this means d in sum(Rd x Jd) (the cartesian product)
#                 # or there is r in Rd and j in Jd such that d = r + j
#                 comb_keys = defaultdict(set)
#                 for r, j in product(Rd_keys, Jd_keys):
#                     comb_keys[r + j].add((r, j))
#                 for d, rj_set in comb_keys.items():
#                     f[R][S][TR | TJ][d] = min(
#                         [
#                             (
#                                 f[R][S][TR][r][0] + f[J][S][TJ][j][0],
#                                 f[R][S][TR][r][1] | f[J][S][TJ][j][1],
#                             )
#                             for r, j in rj_set
#                         ],
#                         key=lambda x: x[0],
#                     )
#
#     stack: tp.List[Bag] = [root]
#     # f is a 4 dimensional matrix:
#     #   1st dimension is Bags (restriction tree vertices)
#     #   2nd dimenstion is possible sets of bps
#     #   3rd dimension is possible sets of variable segments (intervals)
#     #   4th dimension is possible costs
#     f: tp.Dict[
#         Bag,
#         tp.Dict[frozenset, tp.Dict[frozenset, tp.Dict[int, tp.Tuple[int, frozenset]]]],
#     ] = nested_defaultdict(dict, depth=3)
#
#     visited: tp.Set[Bag] = set()
#
#     while stack:
#         R = stack[-1]
#         if bag_type[R] == BagType.Leaf:
#             visited.add(R)
#             stack.pop()
#             leaf(R)
#         else:
#             if R in visited:
#                 # Finished processing R's children its time to process R
#                 # If R is Introduce
#                 if bag_type[R] == BagType.Introduce:
#                     introduce(R)
#                 # If R is forget
#                 elif bag_type[R] == BagType.Forget:
#                     forget(R)
#                     pass
#                 # if R is Join
#                 elif bag_type[R] == BagType.Join:
#                     join(R)
#                     pass
#                 stack.pop()
#                 pass
#             else:
#                 # R is not leaf - push R's children to the stack and mark it
#                 # as visited
#                 stack.extend(list(Rtree[R]))
#                 visited.add(R)
#     return f
#
#
# def MCSCliqueTree(G: nx.Graph):
#     """
#     Ref: Galinier et al. - Chordal graphs and their clique graphs (2005)
#     Data: A graph G = (V, E)
#     Result: If the input graph is chordal: a PEO and an associated
#     clique-tree T = (I, F) where I is the set of maximal cliques
#     begin
#         each vertex of X is initialized with the empty set
#         previousmark = -1
#         j=0
#         for i=n to 1 do
#             choose a vertex x not yet numbered such that |mark(x)| is maximum
#             if mark(x) < previousmark then
#                 j=j+1
#                 create the maximal clique Cj = M(x) U {x}
#                 create the tie between Cj and C(last(x))
#             else
#                 Cj = Cj U {x}
#             for each y neighbour of x do
#                 M(y) = M(y) U {x}
#                 mark(y) = mark(y) + 1
#                 last(y) = x
#             previousmark = mark(x)
#             x is numbered by i
#             C(x) = j
#     end
#     """
#     heap = dict((n, 0) for n in G.nodes)
#     T: tp.List[tp.Tuple[int, int]] = []
#     alpha: tp.Dict[tp.Any, int] = dict()
#     C: tp.Dict[int, tp.Set] = defaultdict(set)
#     C_dict: tp.Dict[tp.Any, int] = dict((n, 0) for n in G.nodes)
#     M: tp.Dict[tp.Any, tp.Set] = dict((n, set()) for n in G.nodes)
#     previousmark = -1
#     j = 0
#     last: tp.Dict[tp.Any, int] = dict((n, 0) for n in G.nodes)
#     for i in range(len(G.nodes), 0, -1):
#         u = max(heap, key=heap.get)
#         marku = heap.pop(u)
#         if marku <= previousmark:
#             j += 1
#             C[j] = M[u] | set([u])
#             T.append((j, C_dict[last[u]]))
#         else:
#             C[j] = C[j] | set([u])
#         for v in G[u]:
#             try:
#                 heap[v] += 1
#                 M[v] = M[v] | set([u])
#                 last[v] = u
#             except KeyError:
#                 pass
#         previousmark = marku
#         alpha[u] = i
#         C_dict[u] = j
#     return T, C, alpha
#
#
# def tree_struct_to_tree(
#     tree_struct: tp.List[tp.Tuple[int, int]], cliques: tp.Dict[int, set]
# ):
#     G = nx.Graph()
#     edges = [(frozenset(cliques[i]), frozenset(cliques[j])) for i, j in tree_struct]
#     G.add_edges_from(edges)
#     G.add_nodes_from(map(lambda x: frozenset(x), cliques.values()))
#     return nx.dfs_tree(G, frozenset(next(iter(cliques.values()))))
#
#
# def algch(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[tp.FrozenSet, tp.Dict[int, tp.Tuple[int, frozenset]]],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     dummy_node: tp.FrozenSet[Gate],
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][0] = (0, frozenset())
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             f[node][v_set][profit(v_set)] = (1, v_set)
#     else:
#         for child in tree[node]:
#             algch(tree, child, f, profit, dummy_node)
#         children = iter(tree[node])
#         child = next(children)
#         node_child_sepset = node.intersection(child)
#         for v in node - dummy_node:
#
#             v_set = frozenset([v])
#             if v in node_child_sepset:
#                 for d, val in f[child][v_set].items():
#                     f[node][v_set][d] = val
#             else:
#                 profits = defaultdict(list)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#                 for i in child_sets:
#                     for weight, gateset in f[child][i].values():
#                         profits[profit(v_set | gateset)].append(
#                             (weight + 1, gateset | v_set)
#                         )
#                 for d in profits.keys():
#                     f[node][v_set][d] = min(profits[d], key=lambda x: x[0])
#         profits = defaultdict(list)
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for d, val in f[child][i].items():
#                 profits[d].append(val)
#         for d in profits.keys():
#             f[node][frozenset()][d] = min(profits[d], key=lambda x: x[0])
#         for child in tqdm.tqdm(children):
#             for v in node - dummy_node:
#                 v_set = frozenset([v])
#                 node_child_sepset = node.intersection(child)
#                 if v in node_child_sepset:
#                     profits = defaultdict(list)
#                     for i, j in product(
#                         f[node][v_set].values(), f[child][v_set].values()
#                     ):
#                         gateset = i[1] | j[1]
#                         profits[profit(gateset)].append((i[0] + j[0] - 1, gateset))
#                     for d in profits.keys():
#                         f[node][v_set][d] = min(profits[d], key=lambda x: x[0])
#                 else:
#                     profits = defaultdict(list)
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for i in child_sets:
#                         for j, k in product(
#                             f[node][v_set].values(), f[child][i].values()
#                         ):
#                             gateset = j[1] | k[1]
#                             profits[profit(gateset)].append((j[0] + k[0], gateset))
#                     for d in profits.keys():
#                         f[node][v_set][d] = min(profits[d], key=lambda x: x[0])
#             profits = defaultdict(list)
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             for i in child_sets:
#                 for j, k in product(
#                     f[node][frozenset()].values(), f[child][i].values()
#                 ):
#                     gateset = j[1] | k[1]
#                     profits[profit(gateset)].append((j[0] + k[0], gateset))
#             for d in profits.keys():
#                 f[node][frozenset()][d] = min(profits[d], key=lambda x: x[0])
#     return f
#
#
# def algch2(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[tp.FrozenSet, tp.Dict[int, tp.Dict[int, tp.Tuple[int, frozenset]]]],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     is_valid_gates: tp.Callable[[Iterable[Gate]], int],
#     dummy_node: tp.FrozenSet[Gate],
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][(0, 0, 0, 0)][0] = (0, frozenset())
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             f[node][v_set][gate_occupancy(v_set)][profit(v_set)] = (1, v_set)
#     else:
#         for child in tree[node]:
#             algch2(tree, child, f, profit, is_valid_gates, dummy_node)
#         children = list(tree[node])
#         child = children[0]
#         node_child_sepset = node.intersection(child)
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             if v in node_child_sepset:
#                 for v_occ in f[child][v_set].keys():
#                     for p_key, p_val in f[child][v_set][v_occ].items():
#                         f[node][v_set][v_occ][p_key] = p_val
#             else:
#                 v_occ = gate_occupancy(v_set)
#                 profits = defaultdict(list)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#                 for i in child_sets:
#                     for occ in f[child][i].keys():
#                         if not is_comp_occ(v_occ, occ):
#                             continue
#                         for weight, gateset in f[child][i][occ].values():
#                             if is_valid_gates(gateset | v_set):
#                                 profits[profit(v_set | gateset)].append(
#                                     (
#                                         weight + 1,
#                                         gateset | v_set,
#                                         gate_occupancy(gateset | v_set),
#                                     )
#                                 )
#                 for d in profits.keys():
#                     profit_min_weight = min(profits[d])
#                     f[node][v_set][profit_min_weight[2]][d] = profit_min_weight[:2]
#         profits = defaultdict(list)
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for occ in f[child][i].keys():
#                 for d, val in f[child][i][occ].items():
#                     profits[d].append(val + (occ,))
#         for d in profits.keys():
#             profit_min_weight = min(profits[d])
#             f[node][frozenset()][profit_min_weight[2]][d] = profit_min_weight[:2]
#         if len(children) < 2:
#             return f
#         for child in tqdm.tqdm(children[1:]):
#             for v in node - dummy_node:
#                 v_set = frozenset([v])
#                 node_child_sepset = node.intersection(child)
#                 if v in node_child_sepset:
#                     profits = defaultdict(list)
#                     for v_occ in f[node][v_set].keys():
#                         for c_occ in f[child][v_set].keys():
#                             if not is_comp_occ(v_occ, c_occ):
#                                 continue
#                             for i, j in product(
#                                 f[node][v_set][v_occ].values(),
#                                 f[child][v_set][c_occ].values(),
#                             ):
#                                 gateset = i[1] | j[1]
#                                 if is_valid_gates(gateset):
#                                     profits[profit(gateset)].append(
#                                         (
#                                             i[0] + j[0] - 1,
#                                             gateset,
#                                             gate_occupancy(gateset),
#                                         )
#                                     )
#                     for d in profits.keys():
#                         profit_min_weight = min(profits[d])
#                         f[node][v_set][profit_min_weight[2]][d] = profit_min_weight[:2]
#                 else:
#                     profits = defaultdict(list)
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for v_occ in f[node][v_set].keys():
#                         for i in child_sets:
#                             for c_occ in f[child][i].keys():
#                                 if not is_comp_occ(v_occ, c_occ):
#                                     continue
#                                 for j, k in product(
#                                     f[node][v_set][v_occ].values(),
#                                     f[child][i][c_occ].values(),
#                                 ):
#                                     gateset = j[1] | k[1]
#                                     if is_valid_gates(gateset):
#                                         profits[profit(gateset)].append(
#                                             (
#                                                 j[0] + k[0],
#                                                 gateset,
#                                                 gate_occupancy(gateset),
#                                             )
#                                         )
#                     for d in profits.keys():
#                         profit_min_weight = min(profits[d])
#                         f[node][v_set][profit_min_weight[2]][d] = profit_min_weight[:2]
#             profits = defaultdict(list)
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             for f_occ in f[node][frozenset()].keys():
#                 for i in child_sets:
#                     for c_occ in f[child][i].keys():
#                         if not is_comp_occ(f_occ, c_occ):
#                             continue
#                         for j, k in product(
#                             f[node][frozenset()][f_occ].values(),
#                             f[child][i][c_occ].values(),
#                         ):
#                             gateset = j[1] | k[1]
#                             if is_valid_gates(gateset):
#                                 profits[profit(gateset)].append(
#                                     (j[0] + k[0], gateset, gate_occupancy(gateset))
#                                 )
#             for d in profits.keys():
#                 profit_min_weight = min(profits[d])
#                 f[node][frozenset()][profit_min_weight[2]][d] = profit_min_weight[:2]
#     return f
#
#
# def algch3(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[tp.FrozenSet, tp.Dict[int, tp.Dict[int, tp.Tuple[int, frozenset]]]],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     is_valid_gates: tp.Callable[[Iterable[Gate]], int],
#     dummy_node: tp.FrozenSet[Gate],
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][(0, 0, 0, 0)][frozenset()] = 0
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             f[node][v_set][gate_occupancy(v_set)][v_set] = 1
#     else:
#         for child in tree[node]:
#             algch3(tree, child, f, profit, is_valid_gates, dummy_node)
#         children = list(tree[node])
#         child = children[0]
#         node_child_sepset = node.intersection(child)
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             if v in node_child_sepset:
#                 for c_occ in f[child][v_set].keys():
#                     for gateset, weight in f[child][v_set][c_occ].items():
#                         f[node][v_set][c_occ][gateset] = weight
#             else:
#                 v_occ = gate_occupancy(v_set)
#                 v_dict = nested_defaultdict(dict, 2)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#                 for i in child_sets:
#                     for occ in f[child][i].keys():
#                         if not is_comp_occ(v_occ, occ):
#                             continue
#                         for gateset, weight in f[child][i][occ].items():
#                             if is_valid_gates(gateset | v_set):
#                                 n_occ = gate_occupancy(gateset | v_set)
#                                 v_dict[n_occ][gateset | v_set] = weight + 1
#                 for n_occ in v_dict.keys():
#                     for gateset, weight in v_dict[n_occ].items():
#                         f[node][v_set][n_occ][gateset] = weight
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for i_occ in f[child][i].keys():
#                 for gateset, weight in f[child][i][i_occ].items():
#                     f[node][frozenset()][i_occ][gateset] = weight
#         # del f[child]
#         if len(children) < 2:
#             return f
#         for child in tqdm.tqdm(children[1:]):
#             for v in node - dummy_node:
#                 v_set = frozenset([v])
#                 node_child_sepset = node.intersection(child)
#                 v_dict = nested_defaultdict(dict, 2)
#                 v_occ_keys = list(f[node][v_set].keys())
#                 if v in node_child_sepset:
#                     for v_occ in v_occ_keys:
#                         for c_occ in f[child][v_set].keys():
#                             if not is_comp_occ(v_occ, c_occ):
#                                 continue
#                             for i, j in product(
#                                 f[node][v_set][v_occ].items(),
#                                 f[child][v_set][c_occ].items(),
#                             ):
#                                 gateset = i[0] | j[0]
#                                 if is_valid_gates(gateset):
#                                     n_occ = gate_occupancy(gateset | v_set)
#
#                 else:
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for v_occ in v_occ_keys:
#                         for i in child_sets:
#                             for c_occ in f[child][i].keys():
#                                 if not is_comp_occ(v_occ, c_occ):
#                                     continue
#                                 for j, k in product(
#                                     f[node][v_set][v_occ].items(),
#                                     f[child][i][c_occ].items(),
#                                 ):
#                                     gateset = j[0] | k[0]
#                                     if is_valid_gates(gateset):
#                                         n_occ = gate_occupancy(gateset | v_set)
#                                         v_dict[n_occ][gateset] = j[1] + k[1]
#                 for n_occ in v_dict.keys():
#                     for gateset, weight in v_dict[n_occ].items():
#                         f[node][v_set][n_occ][gateset] = weight
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             v_keys = list(f[node][frozenset()].keys())
#             v_dict = nested_defaultdict(dict, 2)
#             for v_occ in v_keys:
#                 for i in child_sets:
#                     for c_occ in f[child][i].keys():
#                         if not is_comp_occ(v_occ, c_occ):
#                             continue
#                         for j, k in product(
#                             f[node][frozenset()][v_occ].items(),
#                             f[child][i][c_occ].items(),
#                         ):
#                             gateset = j[0] | k[0]
#                             if is_valid_gates(gateset):
#                                 n_occ = gate_occupancy(gateset)
#                                 v_dict[n_occ][gateset] = j[1] + k[1]
#             for n_occ in v_dict.keys():
#                 for gateset, weight in v_dict[n_occ].items():
#                     f[node][frozenset()][n_occ][gateset] = weight
#     #            del f[child]
#     return f
#
#
# def algch4(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[
#             tp.FrozenSet,
#             tp.Dict[tp.Tuple, tp.Dict[int, tp.Set[tp.Tuple[int, frozenset]]]],
#         ],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     is_valid_gates: tp.Callable[[Iterable[Gate]], int],
#     dummy_node: tp.FrozenSet[Gate],
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#
#     def copy_profits(profits, v_dict):
#         for d in profits.keys():
#             sort_d_profits = sorted(profits[d])
#             profit_min_weight = sort_d_profits[0][0]
#             for sol in sort_d_profits:
#                 if sol[0] > profit_min_weight:
#                     break
#                 v_dict[sol[2]][d].add(sol[:2])
#
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][(0, 0, 0, 0)][0] = set([(0, frozenset())])
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             f[node][v_set][gate_occupancy(v_set)][profit(v_set)] = set([(1, v_set)])
#     else:
#         for child in tree[node]:
#             algch4(tree, child, f, profit, is_valid_gates, dummy_node)
#         children = list(tree[node])
#         child = children[0]
#         node_child_sepset = node.intersection(child)
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             if v in node_child_sepset:
#                 for v_occ in f[child][v_set].keys():
#                     for p_key, p_val in f[child][v_set][v_occ].items():
#                         f[node][v_set][v_occ][p_key] = p_val
#             else:
#                 v_occ = gate_occupancy(v_set)
#                 profits = defaultdict(set)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#
#                 for i in child_sets:
#                     for occ in f[child][i].keys():
#                         if not is_comp_occ(v_occ, occ):
#                             continue
#                         for p_list in f[child][i][occ].values():
#                             for weight, gateset in p_list:
#                                 if is_valid_gates(gateset | v_set):
#                                     profits[profit(v_set | gateset)].add(
#                                         (
#                                             weight + 1,
#                                             gateset | v_set,
#                                             gate_occupancy(gateset | v_set),
#                                         )
#                                     )
#
#                 copy_profits(profits, f[node][v_set])
#
#         profits = defaultdict(set)
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for occ in f[child][i].keys():
#                 for p_list in f[child][i][occ].values():
#                     for d, p_list in f[child][i][occ].items():
#                         for weight, gateset in p_list:
#                             profits[d].add((weight, gateset, occ))
#
#         copy_profits(profits, f[node][frozenset()])
#         del f[child]
#         if len(children) < 2:
#             return f
#
#         for child in tqdm.tqdm(children[1:]):
#             for v in node - dummy_node:
#                 v_set = frozenset([v])
#                 node_child_sepset = node.intersection(child)
#                 profits = defaultdict(set)
#                 if v in node_child_sepset:
#                     for v_occ in f[node][v_set].keys():
#                         for c_occ in f[child][v_set].keys():
#                             if not is_comp_occ(v_occ, c_occ):
#                                 continue
#                             for i, j in product(
#                                 chain.from_iterable(f[node][v_set][v_occ].values()),
#                                 chain.from_iterable(f[child][v_set][c_occ].values()),
#                             ):
#                                 gateset = i[1] | j[1]
#                                 if is_valid_gates(gateset):
#                                     profits[profit(gateset)].add(
#                                         (
#                                             i[0] + j[0] - 1,
#                                             gateset,
#                                             gate_occupancy(gateset),
#                                         )
#                                     )
#
#                 else:
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for v_occ in f[node][v_set].keys():
#                         for i in child_sets:
#                             for c_occ in f[child][i].keys():
#                                 if not is_comp_occ(v_occ, c_occ):
#                                     continue
#                                 for j, k in product(
#                                     chain.from_iterable(f[node][v_set][v_occ].values()),
#                                     chain.from_iterable(f[child][i][c_occ].values()),
#                                 ):
#                                     gateset = j[1] | k[1]
#                                     if is_valid_gates(gateset):
#                                         profits[profit(gateset)].add(
#                                             (
#                                                 j[0] + k[0],
#                                                 gateset,
#                                                 gate_occupancy(gateset),
#                                             )
#                                         )
#
#                 copy_profits(profits, f[node][v_set])
#             profits = defaultdict(set)
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             for v_occ in f[node][frozenset()].keys():
#                 for i in child_sets:
#                     for c_occ in f[child][i].keys():
#                         if not is_comp_occ(v_occ, c_occ):
#                             continue
#                         for j, k in product(
#                             chain.from_iterable(f[node][frozenset()][v_occ].values()),
#                             chain.from_iterable(f[child][i][c_occ].values()),
#                         ):
#                             gateset = j[1] | k[1]
#                             if is_valid_gates(gateset):
#                                 profits[profit(gateset)].add(
#                                     (j[0] + k[0], gateset, gate_occupancy(gateset))
#                                 )
#
#             copy_profits(profits, f[node][frozenset()])
#             del f[child]
#     return f
#
#
# def algch5(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[tp.FrozenSet, tp.Dict[int, tp.Set[tp.Tuple[int, frozenset]]]],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     is_valid_gates: tp.Callable[[Iterable[Gate]], int],
#     dummy_node: tp.FrozenSet[Gate],
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#
#     def copy_profits(profits, v_dict):
#         for d in profits.keys():
#             sort_d_profits = sorted(profits[d])
#             profit_min_weight = sort_d_profits[0][0]
#             for sol in sort_d_profits:
#                 if sol[0] > profit_min_weight:
#                     break
#                 v_dict[d].add(sol)
#
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][0] = set([(0, frozenset())])
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             f[node][v_set][profit(v_set)] = set([(1, v_set)])
#     else:
#         for child in tree[node]:
#             algch5(tree, child, f, profit, is_valid_gates, dummy_node)
#         children = list(tree[node])
#         child = children[0]
#         node_child_sepset = node.intersection(child)
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             if v in node_child_sepset:
#                 for p_key, p_val in f[child][v_set].items():
#                     f[node][v_set][p_key] = p_val
#             else:
#                 profits = defaultdict(set)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#
#                 for i in child_sets:
#                     for p_set in f[child][i].values():
#                         for weight, gateset in p_set:
#                             if is_valid_gates(gateset | v_set):
#                                 profits[profit(v_set | gateset)].add(
#                                     (weight + 1, gateset | v_set)
#                                 )
#
#                 copy_profits(profits, f[node][v_set])
#
#         profits = defaultdict(set)
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for d, p_set in f[child][i].items():
#                 for weight, gateset in p_set:
#                     profits[d].add((weight, gateset))
#
#         copy_profits(profits, f[node][frozenset()])
#         del f[child]
#         if len(children) < 2:
#             return f
#
#         for child in tqdm.tqdm(children[1:]):
#             for v in node - dummy_node:
#                 v_set = frozenset([v])
#                 node_child_sepset = node.intersection(child)
#                 profits = defaultdict(set)
#                 if v in node_child_sepset:
#                     for i, j in product(
#                         chain.from_iterable(f[node][v_set].values()),
#                         chain.from_iterable(f[child][v_set].values()),
#                     ):
#                         gateset = i[1] | j[1]
#                         if is_valid_gates(gateset):
#                             profits[profit(gateset)].add((i[0] + j[0] - 1, gateset))
#
#                 else:
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for j, k in product(
#                         chain.from_iterable(f[node][v_set].values()),
#                         chain.from_iterable(f[child][i].values()),
#                     ):
#                         gateset = j[1] | k[1]
#                         if is_valid_gates(gateset):
#                             profits[profit(gateset)].add((j[0] + k[0], gateset))
#
#                 copy_profits(profits, f[node][v_set])
#             profits = defaultdict(set)
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             for j, k in product(
#                 chain.from_iterable(f[node][frozenset()].values()),
#                 chain.from_iterable(f[child][i].values()),
#             ):
#                 gateset = j[1] | k[1]
#                 if is_valid_gates(gateset):
#                     profits[profit(gateset)].add((j[0] + k[0], gateset))
#
#             copy_profits(profits, f[node][frozenset()])
#             del f[child]
#     return f
#
#
# def algch6(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[tp.FrozenSet, tp.Dict[int, tp.Set[tp.Tuple[int, frozenset]]]],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     is_valid_gates: tp.Callable[[Iterable[Gate]], int],
#     dummy_node: tp.FrozenSet[Gate],
#     C: int = 12,
#     empty_cost: int = 618624,
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#
#     def copy_profits(profits, v_dict):
#         for d in profits.keys():
#             v_dict[d] = v_dict[d] | profits[d]
#             # sort_d_profits = sorted(profits[d])
#             # profit_min_weight = sort_d_profits[0][0]
#             # for sol in sort_d_profits:
#             #     if sol[0] > profit_min_weight:
#             #         break
#             #     v_dict[profit_min_weight].add(sol)
#
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][0] = set([(0, empty_cost, frozenset())])
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             f[node][v_set][1] = set([(1, profit(v_set), v_set)])
#     else:
#         for child in tree[node]:
#             algch6(tree, child, f, profit, is_valid_gates, dummy_node)
#         children = list(tree[node])
#         child = children[0]
#         node_child_sepset = node.intersection(child)
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             if v in node_child_sepset:
#                 copy_profits(f[child][v_set], f[node][v_set])
#             else:
#                 profits = defaultdict(set)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#
#                 for i in child_sets:
#                     for weight, p_set in f[child][i].items():
#                         if weight + 1 > C:
#                             continue
#                         for _, p, gateset in p_set:
#                             if (
#                                 is_valid_gates(gateset | v_set)
#                                 and profit(gateset | v_set) <= p
#                             ):
#                                 profits[weight + 1].add(
#                                     (
#                                         weight + 1,
#                                         profit(gateset | v_set),
#                                         gateset | v_set,
#                                     )
#                                 )
#                 copy_profits(profits, f[node][v_set])
#         profits = defaultdict(set)
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for p, p_set in f[child][i].items():
#                 profits[p] |= p_set
#
#         copy_profits(profits, f[node][frozenset()])
#
#         del f[child]
#         if len(children) < 2:
#             return f
#
#         for child in tqdm.tqdm(children[1:]):
#             for v in node - dummy_node:
#                 v_set = frozenset([v])
#                 node_child_sepset = node.intersection(child)
#                 profits = defaultdict(set)
#                 if v in node_child_sepset:
#                     for w_i, w_j in product(
#                         f[node][v_set].keys(), f[child][v_set].keys()
#                     ):
#                         if w_i + w_j > C:
#                             continue
#                         for i, j in product(f[node][v_set][w_i], f[child][v_set][w_j]):
#                             gateset = i[-1] | j[-1]
#                             if (
#                                 is_valid_gates(gateset)
#                                 and profit(gateset) <= i[1]
#                                 and profit(gateset) <= j[1]
#                             ):
#                                 profits[len(gateset)].add(
#                                     (len(gateset), profit(gateset), gateset)
#                                 )
#
#                 else:
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for i in child_sets:
#                         for w_i, w_j in product(
#                             f[node][v_set].keys(), f[child][i].keys()
#                         ):
#                             if w_i + w_j > C:
#                                 continue
#                             for j, k in product(f[node][v_set][w_i], f[child][i][w_j]):
#                                 gateset = j[-1] | k[-1]
#                                 if (
#                                     is_valid_gates(gateset)
#                                     and profit(gateset) <= j[1]
#                                     and profit(gateset) <= k[1]
#                                 ):
#                                     profits[len(gateset)].add(
#                                         (len(gateset), profit(gateset), gateset)
#                                     )
#
#                 copy_profits(profits, f[node][v_set])
#
#             profits = defaultdict(set)
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             for i in child_sets:
#                 for w_i, w_j in product(
#                     f[node][frozenset()].keys(), f[child][i].keys()
#                 ):
#                     if w_i + w_j > C:
#                         continue
#                     for j, k in product(f[node][frozenset()][w_i], f[child][i][w_j]):
#                         gateset = j[-1] | k[-1]
#                         if (
#                             is_valid_gates(gateset)
#                             and profit(gateset) <= j[1]
#                             and profit(gateset) <= k[1]
#                         ):
#                             profits[len(gateset)].add(
#                                 (len(gateset), profit(gateset), gateset)
#                             )
#
#             copy_profits(profits, f[node][frozenset()])
#             del f[child]
#     return f
#
#
# def algch7(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[tp.FrozenSet, tp.Dict[int, tp.Set[tp.Tuple[int, int, frozenset]]]],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     is_valid_gates: tp.Callable[[Iterable[Gate]], int],
#     dummy_node: tp.FrozenSet[Gate],
#     C: int = 12,
#     empty_cost: int = 618624,
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#
#     def copy_profits(profits, v_dict):
#         for d in profits.keys():
#             v_dict[d] = v_dict[d] | profits[d]
#             # sort_d_profits = sorted(profits[d])
#             # try:
#             #     profit_min_weight = sort_d_profits[0][0]
#             #     for sol in sort_d_profits:
#             #         if sol[0] > profit_min_weight:
#             #             break
#             #         v_dict[profit_min_weight].add(sol)
#             # except IndexError:
#             #     pass
#
#     def remove_unnecessary(f_dict, node):
#         d_max: tp.Dict[int, int] = {}
#         for v_set in f_dict[node].keys():
#             for d, d_set in f_dict[node][v_set].items():
#                 try:
#                     d_max[d] = max([ss[1] for ss in d_set])
#                 except ValueError:
#                     d_max[d] = np.inf
#         for v_set in f_dict[node].keys():
#             v_set_dict = defaultdict(set)
#             for d, d_set in f_dict[node][v_set].items():
#                 try:
#                     v_set_dict[d] = set(filter(lambda ss: ss[1] <= d_max[d - 1], d_set))
#                 except KeyError:
#                     v_set_dict[d] = d_set
#             f_dict[node][v_set] = v_set_dict
#
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][0] = set([(0, empty_cost, frozenset())])
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             f[node][v_set][1] = set([(1, profit(v_set), v_set)])
#     else:
#         for child in tree[node]:
#             algch7(tree, child, f, profit, is_valid_gates, dummy_node)
#         children = list(tree[node])
#         child = children[0]
#         node_child_sepset = node.intersection(child)
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             if v in node_child_sepset:
#                 copy_profits(f[child][v_set], f[node][v_set])
#             else:
#                 profits = defaultdict(set)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#
#                 for i in child_sets:
#                     for weight, p_set in f[child][i].items():
#                         if weight + 1 > C:
#                             continue
#                         for _, p, gateset in p_set:
#                             if (
#                                 is_valid_gates(gateset | v_set)
#                                 and profit(gateset | v_set) <= p
#                             ):
#                                 profits[weight + 1].add(
#                                     (
#                                         weight + 1,
#                                         profit(gateset | v_set),
#                                         gateset | v_set,
#                                     )
#                                 )
#                 copy_profits(profits, f[node][v_set])
#         profits = defaultdict(set)
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for p, p_set in f[child][i].items():
#                 profits[p] |= p_set
#
#         copy_profits(profits, f[node][frozenset()])
#
#         del f[child]
#         remove_unnecessary(f, node)
#         if len(children) < 2:
#             return f
#
#         for child in tqdm.tqdm(children[1:]):
#             for v in node - dummy_node:
#                 v_set = frozenset([v])
#                 node_child_sepset = node.intersection(child)
#                 profits = defaultdict(set)
#                 if v in node_child_sepset:
#                     for w_i, w_j in product(
#                         f[node][v_set].keys(), f[child][v_set].keys()
#                     ):
#                         if w_i + w_j > C:
#                             continue
#                         for i, j in product(f[node][v_set][w_i], f[child][v_set][w_j]):
#                             gateset = i[-1] | j[-1]
#                             if (
#                                 is_valid_gates(gateset)
#                                 and profit(gateset) <= i[1]
#                                 and profit(gateset) <= j[1]
#                             ):
#                                 profits[len(gateset)].add(
#                                     (len(gateset), profit(gateset), gateset)
#                                 )
#
#                 else:
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for i in child_sets:
#                         for w_i, w_j in product(
#                             f[node][v_set].keys(), f[child][i].keys()
#                         ):
#                             if w_i + w_j > C:
#                                 continue
#                             for j, k in product(f[node][v_set][w_i], f[child][i][w_j]):
#                                 gateset = j[-1] | k[-1]
#                                 if (
#                                     is_valid_gates(gateset)
#                                     and profit(gateset) <= j[1]
#                                     and profit(gateset) <= k[1]
#                                 ):
#                                     profits[len(gateset)].add(
#                                         (len(gateset), profit(gateset), gateset)
#                                     )
#
#                 copy_profits(profits, f[node][v_set])
#
#             profits = defaultdict(set)
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             for i in child_sets:
#                 for w_i, w_j in product(
#                     f[node][frozenset()].keys(), f[child][i].keys()
#                 ):
#                     if w_i + w_j > C:
#                         continue
#                     for j, k in product(f[node][frozenset()][w_i], f[child][i][w_j]):
#                         gateset = j[-1] | k[-1]
#                         if (
#                             is_valid_gates(gateset)
#                             and profit(gateset) <= j[1]
#                             and profit(gateset) <= k[1]
#                         ):
#                             profits[len(gateset)].add(
#                                 (len(gateset), profit(gateset), gateset)
#                             )
#
#             copy_profits(profits, f[node][frozenset()])
#             del f[child]
#             remove_unnecessary(f, node)
#     return f
#
#
# def algch8(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[
#             tp.FrozenSet,
#             tp.Dict[int, tp.Dict[int, tp.Set[tp.Tuple[int, int, frozenset]]]],
#         ],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     is_valid_gates: tp.Callable[[Iterable[Gate]], int],
#     dummy_node: tp.FrozenSet[Gate],
#     C: int = 6,
#     max_P: int = 2143,
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#
#     def copy_profits(profits, v_dict):
#         for d in profits.keys():
#             # v_dict[d] = frozenset([min(profits[d])])
#             sort_d_profits = sorted(profits[d])
#             try:
#                 profit_min_weight = sort_d_profits[0][0]
#                 for sol in sort_d_profits:
#                     if sol[0] > profit_min_weight:
#                         break
#                     v_dict[d].add(sol)
#             except IndexError:
#                 pass
#
#     def copy_node(child, node):
#         for p, p_set in child.items():
#             node[p] = p_set
#
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][0] = set([(0, frozenset())])
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             f[node][v_set][profit(v_set)] = set([(1, v_set)])
#
#     else:
#         for child in tree[node]:
#             algch8(tree, child, f, profit, is_valid_gates, dummy_node)
#         children = list(tree[node])
#         child = children[0]
#         node_child_sepset = node.intersection(child)
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             v_prof = profit(v_set)
#             if v in node_child_sepset:
#                 copy_node(f[child][v_set], f[node][v_set])
#             else:
#                 profits = defaultdict(set)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#
#                 for i in child_sets:
#                     for p, p_set in f[child][i].items():
#                         # if profit is more than max or the weight of the the
#                         # current set does not permit addition - continue
#                         if p + v_prof > max_P or next(iter(p_set))[0] + 1 > C:
#                             # print('skipping solution')
#                             continue
#                         for weight, gateset in p_set:
#                             gateset = gateset | v_set
#                             if is_valid_gates(gateset):
#                                 profits[p + v_prof].add((weight + 1, gateset | v_set))
#                 copy_profits(profits, f[node][v_set])
#
#         profits = defaultdict(set)
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for p, p_set in f[child][i].items():
#                 profits[p] |= p_set
#
#         copy_profits(profits, f[node][frozenset()])
#
#         if len(children) < 2:
#             return f
#
#         for child in tqdm.tqdm(children[1:]):
#             for v in node - dummy_node:
#                 v_set = frozenset([v])
#                 v_prof = profit(v_set)
#                 node_child_sepset = node.intersection(child)
#                 profits = defaultdict(set)
#                 if v in node_child_sepset:
#                     for p_j, p_k in product(
#                         f[node][v_set].keys(), f[child][v_set].keys()
#                     ):
#                         if (
#                             p_j + p_k > max_P
#                             or next(iter(f[node][v_set][p_j]))[0]
#                             + next(iter(f[child][v_set][p_k]))[0]
#                             > C
#                         ):
#                             # print('skipping solution')
#                             continue
#                         for j, k in product(
#                             (f[node][v_set][p_j]), (f[child][v_set][p_k])
#                         ):
#
#                             gateset = j[1] | k[1]
#                             if is_valid_gates(gateset):
#                                 profits[p_j + p_k - v_prof].add(
#                                     (j[0] + k[0] - 1, gateset)
#                                 )
#
#                 else:
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for i in child_sets:
#                         for p_j, p_k in product(
#                             f[node][v_set].keys(), f[child][i].keys()
#                         ):
#                             if (
#                                 p_j + p_k > max_P
#                                 or next(iter(f[node][v_set][p_j]))[0]
#                                 + next(iter(f[child][i][p_k]))[0]
#                                 > C
#                             ):
#                                 # print('skipping solution')
#                                 continue
#                             for j, k in product(
#                                 (f[node][v_set][p_j]), (f[child][i][p_k])
#                             ):
#                                 gateset = j[1] | k[1]
#                                 if is_valid_gates(gateset):
#                                     profits[p_j + p_k].add((j[0] + k[0], gateset))
#
#                 copy_profits(profits, f[node][v_set])
#
#             profits = defaultdict(set)
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             v_set = frozenset()
#             for i in child_sets:
#                 for p_j, p_k in product(f[node][v_set].keys(), f[child][i].keys()):
#                     if (
#                         p_j + p_k > max_P
#                         or next(iter(f[node][v_set][p_j]))[0]
#                         + next(iter(f[child][i][p_k]))[0]
#                         > C
#                     ):
#                         # print('skipping solution')
#                         continue
#                     for j, k in product((f[node][v_set][p_j]), (f[child][i][p_k])):
#                         gateset = j[1] | k[1]
#                         if is_valid_gates(gateset):
#                             profits[p_j + p_k].add((j[0] + k[0], gateset))
#
#             copy_profits(profits, f[node][v_set])
#     return f
#
#
# def algch9(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[
#             tp.FrozenSet,
#             tp.Dict[int, tp.Dict[int, tp.Set[tp.Tuple[int, int, frozenset]]]],
#         ],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     C: int = 12,
#     max_P: int = 2143,
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#
#     def copy_profits(profits, v_dict):
#         for d in profits.keys():
#             v_dict[d] = frozenset([min(profits[d])])
#         #     sort_d_profits = sorted(profits[d])
#         #     try:
#         #         profit_min_weight = sort_d_profits[0][0]
#         #         for sol in sort_d_profits:
#         #             if sol[0] > profit_min_weight:
#         #                 break
#         #             v_dict[d].add(sol)
#         #     except IndexError:
#         #         pass
#
#     def copy_node(child, node):
#         for p, p_set in child.items():
#             node[p] = p_set
#
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][0] = set([(0, frozenset())])
#         for v in node:
#             v_set = frozenset([v])
#             f[node][v_set][profit(v_set)] = set([(1, v_set)])
#
#     else:
#         for child in tree[node]:
#             algch9(tree, child, f, profit)
#         children = list(tree[node])
#         child = children[0]
#         node_child_sepset = node.intersection(child)
#         for v in node:
#             v_set = frozenset([v])
#             v_prof = profit(v_set)
#             if v in node_child_sepset:
#                 copy_node(f[child][v_set], f[node][v_set])
#             else:
#                 profits = defaultdict(set)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#
#                 for i in child_sets:
#                     for p, p_set in f[child][i].items():
#                         # if profit is more than max or the weight of the the
#                         # current set does not permit addition - continue
#                         if p + v_prof > max_P or next(iter(p_set))[0] + 1 > C:
#                             continue
#                         for weight, gateset in p_set:
#                             gateset = gateset | v_set
#                             profits[p + v_prof].add((weight + 1, gateset))
#                 copy_profits(profits, f[node][v_set])
#
#         profits = defaultdict(set)
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for p, p_set in f[child][i].items():
#                 profits[p] |= p_set
#
#         copy_profits(profits, f[node][frozenset()])
#
#         del f[child]
#         if len(children) < 2:
#             return f
#
#         for child in tqdm.tqdm(children[1:]):
#             for v in node:
#                 v_set = frozenset([v])
#                 v_prof = profit(v_set)
#                 node_child_sepset = node.intersection(child)
#                 profits = defaultdict(set)
#                 if v in node_child_sepset:
#                     for p_j, p_k in product(
#                         f[node][v_set].keys(), f[child][v_set].keys()
#                     ):
#                         if (
#                             p_j + p_k > max_P
#                             or next(iter(f[node][v_set][p_j]))[0]
#                             + next(iter(f[child][v_set][p_k]))[0]
#                             > C
#                         ):
#                             continue
#                         for j, k in product(
#                             (f[node][v_set][p_j]), (f[child][v_set][p_k])
#                         ):
#
#                             gateset = j[1] | k[1]
#                             profits[p_j + p_k - v_prof].add((j[0] + k[0] - 1, gateset))
#
#                 else:
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for i in child_sets:
#                         for p_j, p_k in product(
#                             f[node][v_set].keys(), f[child][i].keys()
#                         ):
#                             if (
#                                 p_j + p_k > max_P
#                                 or next(iter(f[node][v_set][p_j]))[0]
#                                 + next(iter(f[child][i][p_k]))[0]
#                                 > C
#                             ):
#                                 continue
#                             for j, k in product(
#                                 (f[node][v_set][p_j]), (f[child][i][p_k])
#                             ):
#                                 gateset = j[1] | k[1]
#                                 profits[p_j + p_k].add((j[0] + k[0], gateset))
#
#                 copy_profits(profits, f[node][v_set])
#
#             profits = defaultdict(set)
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             v_set = frozenset()
#             for i in child_sets:
#                 for p_j, p_k in product(f[node][v_set].keys(), f[child][i].keys()):
#                     if (
#                         p_j + p_k > max_P
#                         or next(iter(f[node][v_set][p_j]))[0]
#                         + next(iter(f[child][i][p_k]))[0]
#                         > C
#                     ):
#                         continue
#                     for j, k in product((f[node][v_set][p_j]), (f[child][i][p_k])):
#                         gateset = j[1] | k[1]
#                         profits[p_j + p_k].add((j[0] + k[0], gateset))
#
#             copy_profits(profits, f[node][v_set])
#             del f[child]
#     return f
#
#
# def create_connecet_edges(G: nx.Graph, dummy_node):
#     connect_edges = []
#     for cc in nx.connected_components(G):
#         min_deg = np.inf
#         min_node = None
#         for node in cc:
#             if G.degree(node) < min_deg:
#                 min_deg = G.degree(node)
#                 min_node = node
#         connect_edges.append((dummy_node, min_node))
#
#     return connect_edges
#
#
# def algch8(
#     tree: nx.DiGraph,
#     node: tp.FrozenSet[Gate],
#     f: tp.Dict[
#         tp.FrozenSet[Gate],
#         tp.Dict[
#             tp.FrozenSet,
#             tp.Dict[int, tp.Dict[int, tp.Set[tp.Tuple[int, int, frozenset]]]],
#         ],
#     ],
#     profit: tp.Callable[[Iterable[Gate]], int],
#     is_valid_gates: tp.Callable[[Iterable[Gate]], int],
#     dummy_node: tp.FrozenSet[Gate],
#     C: int = 6,
#     max_P: int = 2143,
# ):
#     """
#     Ref: Pferschy and Schauser - The knapsack problem with conflict graphs (2009)
#     """
#
#     def copy_profits(profits, v_dict):
#         for d in profits.keys():
#             # v_dict[d] = frozenset([min(profits[d])])
#             sort_d_profits = sorted(profits[d])
#             try:
#                 profit_min_weight = sort_d_profits[0][0]
#                 for sol in sort_d_profits:
#                     if sol[0] > profit_min_weight:
#                         break
#                     v_dict[d].add(sol)
#             except IndexError:
#                 pass
#
#     def copy_node(child, node):
#         for p, p_set in child.items():
#             node[p] = p_set
#
#     if tree.out_degree(node) == 0:
#         f[node][frozenset()][0] = set([(0, frozenset())])
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             f[node][v_set][profit(v_set)] = set([(1, v_set)])
#
#     else:
#         for child in tree[node]:
#             algch8(tree, child, f, profit, is_valid_gates, dummy_node)
#         children = list(tree[node])
#         child = children[0]
#         node_child_sepset = node.intersection(child)
#         for v in node - dummy_node:
#             v_set = frozenset([v])
#             v_prof = profit(v_set)
#             if v in node_child_sepset:
#                 copy_node(f[child][v_set], f[node][v_set])
#             else:
#                 profits = defaultdict(set)
#                 child_sets = list(
#                     map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#                 ) + [frozenset()]
#
#                 for i in child_sets:
#                     for p, p_set in f[child][i].items():
#                         # if profit is more than max or the weight of the the
#                         # current set does not permit addition - continue
#                         if p + v_prof > max_P or next(iter(p_set))[0] + 1 > C:
#                             # print('skipping solution')
#                             continue
#                         for weight, gateset in p_set:
#                             gateset = gateset | v_set
#                             if is_valid_gates(gateset):
#                                 profits[p + v_prof].add((weight + 1, gateset | v_set))
#                 copy_profits(profits, f[node][v_set])
#
#         profits = defaultdict(set)
#         child_sets = list(
#             map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#         ) + [frozenset()]
#         for i in child_sets:
#             for p, p_set in f[child][i].items():
#                 profits[p] |= p_set
#
#         copy_profits(profits, f[node][frozenset()])
#
#         if len(children) < 2:
#             return f
#
#         for child in tqdm.tqdm(children[1:]):
#             for v in node - dummy_node:
#                 v_set = frozenset([v])
#                 v_prof = profit(v_set)
#                 node_child_sepset = node.intersection(child)
#                 profits = defaultdict(set)
#                 if v in node_child_sepset:
#                     for p_j, p_k in product(
#                         f[node][v_set].keys(), f[child][v_set].keys()
#                     ):
#                         if (
#                             p_j + p_k > max_P
#                             or next(iter(f[node][v_set][p_j]))[0]
#                             + next(iter(f[child][v_set][p_k]))[0]
#                             > C
#                         ):
#                             # print('skipping solution')
#                             continue
#                         for j, k in product(
#                             (f[node][v_set][p_j]), (f[child][v_set][p_k])
#                         ):
#
#                             gateset = j[1] | k[1]
#                             if is_valid_gates(gateset):
#                                 profits[p_j + p_k - v_prof].add(
#                                     (j[0] + k[0] - 1, gateset)
#                                 )
#
#                 else:
#                     child_sets = list(
#                         map(
#                             lambda x: frozenset([x]),
#                             child.difference(node_child_sepset),
#                         )
#                     ) + [frozenset()]
#                     for i in child_sets:
#                         for p_j, p_k in product(
#                             f[node][v_set].keys(), f[child][i].keys()
#                         ):
#                             if (
#                                 p_j + p_k > max_P
#                                 or next(iter(f[node][v_set][p_j]))[0]
#                                 + next(iter(f[child][i][p_k]))[0]
#                                 > C
#                             ):
#                                 # print('skipping solution')
#                                 continue
#                             for j, k in product(
#                                 (f[node][v_set][p_j]), (f[child][i][p_k])
#                             ):
#                                 gateset = j[1] | k[1]
#                                 if is_valid_gates(gateset):
#                                     profits[p_j + p_k].add((j[0] + k[0], gateset))
#
#                 copy_profits(profits, f[node][v_set])
#
#             profits = defaultdict(set)
#             child_sets = list(
#                 map(lambda x: frozenset([x]), child.difference(node_child_sepset))
#             ) + [frozenset()]
#             v_set = frozenset()
#             for i in child_sets:
#                 for p_j, p_k in product(f[node][v_set].keys(), f[child][i].keys()):
#                     if (
#                         p_j + p_k > max_P
#                         or next(iter(f[node][v_set][p_j]))[0]
#                         + next(iter(f[child][i][p_k]))[0]
#                         > C
#                     ):
#                         # print('skipping solution')
#                         continue
#                     for j, k in product((f[node][v_set][p_j]), (f[child][i][p_k])):
#                         gateset = j[1] | k[1]
#                         if is_valid_gates(gateset):
#                             profits[p_j + p_k].add((j[0] + k[0], gateset))
#
#             copy_profits(profits, f[node][v_set])
#     return f
#
#
# def LexBFSOrdering(G: nx.Graph, initial: int = 0):
#     alpha = {node: 0 for node in G}
#     L = [set(G.nodes)]
#     i = initial
#     while L:
#         try:
#             x = max(L[-1], key=lambda node: G.degree(node))
#             L[-1].remove(x)
#             alpha[x] = i
#             i -= 1
#             x_neighbours = list(G[x])
#             new_L = []
#             for j, S in enumerate(L):
#                 x_S_neighbours = S.intersection(x_neighbours)
#                 if len(x_S_neighbours) != 0 and x_S_neighbours != S:
#                     new_L.append(S - x_S_neighbours)
#                     new_L.append(x_S_neighbours)
#                 else:
#                     new_L.append(S)
#             L = new_L
#         except ValueError:
#             L.pop()
#
#     return alpha
#
#
# def intersection_graph(G: nx.Graph, elim_order, alpha):
#     tree = nx.Graph()
#     for i, node in enumerate(elim_order):
#         if i + 1 == len(elim_order):
#             continue
#         n_plus = set(elim_order[i + 1 :]).intersection(G[node])
#         parent = min(n_plus, key=alpha.get)
#         tree.add_edge(node, parent)
#     return tree
#
#
# def edge_weight(
#     edge: tp.Tuple[Gate, Gate], var_pos_weight: tp.Dict[tp.Tuple[int, int], float]
# ):
#     acc = 0.0
#     for (seg_start, seg_end), cost in var_pos_weight.items():
#         if edge[0].index < seg_start and seg_end < edge[1].index:
#             acc += cost
#     segment_len = edge[1].index - edge[0].index
#     if acc > 0:
#         if segment_len < 20 or segment_len > 40:
#             acc = np.inf
#     else:
#         if segment_len < 20:
#             acc = np.inf
#     return acc
#
#
# def count_edge_weight(
#     edge: tp.Tuple[Gate, Gate], bps_count: tp.Dict[str, int], incomp_graph: nx.Graph
# ) -> float:
#     return bps_count.get(edge[0].bps, 1) ** (
#         len(incomp_graph[edge[0].bps]) + 1
#     ) * bps_count.get(edge[1].bps, 1) ** (len(incomp_graph[edge[1].bps]) + 1)
#
#
# def prob_edge_weight(
#     edge: tp.Tuple[Gate, Gate],
#     bps_to_edges: tp.Dict[str, tp.Set[Gate]],
#     incomp_graph: nx.Graph,
# ) -> float:
#     incomp_edges = (
#         sum(len(bps_to_edges[gate.bps]) for gate in edge) - 2
#     )  # Remove our given edge
#     incomp_edges += sum(
#         len(bps_to_edges[bps]) for gate in edge for bps in incomp_graph[gate.bps]
#     )
#     total_edges = sum(len(edges) for edges in bps_to_edges.values())
#     return incomp_edges / total_edges
#
#
# def valid_edge(
#     edge: tp.Tuple[Gate, Gate], incomp_graph: nx.Graph, var_pos_arr: np.array
# ) -> bool:
#     if incomp_graph.has_edge(edge[0].bps, edge[1].bps):
#         return False
#     across_var_seg = np.any(
#         np.logical_and(edge[0].idx < var_pos_arr, var_pos_arr < edge[1].idx)
#     )
#     segment_len = edge[1].idx - edge[0].idx
#     if across_var_seg:
#         if segment_len < 20 or segment_len > 36:
#             return False
#     else:
#         if segment_len < 20:
#             return False
#     return True
#
#
# def is_var_edge(edge: tp.Tuple[Gate, Gate], var_pos_arr: np.array) -> bool:
#     return np.any(np.logical_and(edge[0].idx < var_pos_arr, var_pos_arr < edge[1].idx))
#
#
# def nplusi(graph, alpha, elim_order, node):
#     return set(elim_order[alpha[node] :]).intersection(graph[node])
#
#
# def create_default_weight_func(
#     dna_pos_n_codons: Dict[int, List[str]]
# ) -> Callable[[Gate, Gate], int]:
#     var_pos_arr = np.array(list(dna_pos_n_codons.keys()))
#     var_pos_codons = np.array([len(codons) for codons in dna_pos_n_codons.values()])
#
#     @lru_cache(maxsize=67108864)
#     def edge_weight(nd1: Gate, nd2: Gate) -> int:
#         # base_cost = 0
#         # if np.any(np.logical_and(nd1.idx < var_pos_arr, var_pos_arr < nd2.idx)):
#         #     base_cost = 1
#         # return (nd2.idx - nd1.idx) * np.product(
#         #     [
#         #         len(codons) if nd1.idx < pos < nd2.idx else base_cost
#         #         for pos, codons in dna_pos_n_codons.items()
#         #     ]
#         # )
#         return (nd2.idx - nd1.idx) * np.product(
#             var_pos_codons[np.logical_and(nd1.idx < var_pos_arr, var_pos_arr < nd2.idx)]
#         )
#
#     return edge_weight
#
#
# def make_profit_func(
#     dna_pos_n_codons: Dict[int, List[str]], source: Gate, target: Gate
# ) -> Callable[[tp.Iterable[Gate]], int]:
#
#     two_gate_profit = create_default_weight_func(dna_pos_n_codons)
#
#     def profit_func(gates: tp.Iterable[Gate]) -> int:
#         sorted_gates = sorted(gates)
#         ends_profits = 0
#         try:
#             ends_profits += sum(
#                 [
#                     two_gate_profit(source, sorted_gates[0]),
#                     two_gate_profit(sorted_gates[-1], target),
#                 ]
#             )
#         except IndexError:
#             pass
#         profits = (
#             two_gate_profit(g1, g2)
#             for g1, g2 in zip(sorted_gates[:-1], sorted_gates[1:])
#         )
#         return sum(profits) + ends_profits
#
#     return profit_func
#
#
# def enum_is_iter(
#     tree: nx.DiGraph, pre_order: tp.List, i: int, S: set, dummy_node
# ) -> tp.Generator[set, None, None]:
#     stack = [(0, S)]
#     while stack:
#         i, S = stack.pop()
#         if i == len(pre_order):
#             yield S
#         else:
#             if len(pre_order[i] & S) == 0:
#                 parent_v = reduce(
#                     set.union,
#                     [parent for parent in tree.predecessors(pre_order[i])],
#                     set(),
#                 )
#                 for v in pre_order[i] - dummy_node - parent_v:
#                     stack.append((i + 1, S | set([v])))
#             stack.append((i + 1, S))
#
#
# def chordal_mwis(G: nx.Graph, alpha: tp.List, weight_dict: tp.Dict) -> tp.Set:
#     reds: tp.List = []
#     blues: tp.Set = set()
#     for i, v in enumerate(alpha):
#         if weight_dict[v] > 0:
#             reds.append(v)
#             for u in G.subgraph(alpha[i:])[v]:
#                 weight_dict[u] = weight_dict[u] - weight_dict[v]
#                 if weight_dict[u] <= 0:
#                     weight_dict[u] = 0
#             weight_dict[v] = 0
#     for v in reversed(reds):
#         if len(blues.intersection(G[v])) == 0:
#             blues.add(v)
#     return blues
#
#
# def min_end(G: nx.Graph) -> int:
#     min_end = np.NINF
#     for node: tp.Tuple[Gate, Gate] in G:
#         if node[1].idx < min_end:
#             min_end = node[1].idx
#     return min_end
#
#
# def interval_colored_search_tree(
#     G: nx.Graph,
#     k: int,
#     coloring: tp.Dict[tp.Any, tp.Set[int]],
# ) -> tp.Set[tp.Any]:
#     if k == 0:
#         return set()
#     if len(G) == 0:
#         return set()
#     first_end = min_end(G)
#
#
# ggdata = GGData()
# nodes = ggdata.filter_self_binding_gates()
# edges = [
#     (v1, v2)
#     for v1, v2 in combinations(nodes, 2)
#     if ggdata.gates_all_scores(v1, v2) >= 1000
# ]
#
# R_graph = nx.Graph()
# R_graph.add_nodes_from(nodes)
# R_graph.add_edges_from(edges)
# del edges, nodes
#
# chordal_rgraph = nx.complete_to_chordal_graph(R_graph)[0]
#
# dna = parse_dna(
#     "/home/labs/fleishman/shayho/Code/dawdlib/dawdlib/golden_gate/tests/input_921/921_DNA_seq"
# )
# deg_table = pd.read_csv(
#     "/home/labs/fleishman/shayho/Downloads/921.csv",
#     na_filter=True,
#     keep_default_na=False,
# )
# var_poss = expand_dna_var_poss(deg_table[TableColNames.DNA_POS.value].tolist())
# var_pos_arr = np.array(var_poss)
#
# var_segemets = []
# prev_idx = np.NINF
# for var_idx in var_pos_arr:
#     if prev_idx + 4 < var_idx:
#         var_segemets.append([var_idx])
#     else:
#         var_segemets[-1].append(var_idx)
#     prev_idx = var_idx
# var_segemets = pd.DataFrame(var_segemets).values
# var_segemets[np.isnan(var_segemets)] = np.inf
#
# d_graph, src, target = make_default_graph(
#     GraphMaker(ggdata), dna, var_poss, deg_table_to_dict(deg_table), 3, 40, 3
# )
#
# d_graph.remove_nodes_from([src, target])
# var_edges = [edge for edge in d_graph.edges if is_var_edge(edge, var_pos_arr)]
# var_nodes = list(reduce(set.union, [set([edge[0], edge[1]]) for edge in var_edges]))
#
# dr_graph = nx.Graph()
# dr_edges = [
#     (n1, n2)
#     for n1, n2 in combinations(var_nodes, 2)
#     if n1.bps == n2.bps or not is_independent_set(R_graph, set([n1.bps, n2.bps]))
# ]
# dr_graph.add_nodes_from(var_nodes)
# dr_graph.add_edges_from(dr_edges)
#
#
# if not nx.is_chordal(dr_graph):
#     dr_graph = nx.Graph()
#     dr_edges = [
#         (n1, n2)
#         for n1, n2 in combinations(var_nodes, 2)
#         if n1.bps == n2.bps
#         or not is_independent_set(chordal_rgraph, set([n1.bps, n2.bps]))
#     ]
#     dr_graph.add_nodes_from(var_nodes)
#     dr_graph.add_edges_from(dr_edges)
#
# if not nx.is_connected(dr_graph):
#     DUMMY = Gate(np.inf, "DUMMY")
#     connect_edges = create_connecet_edges(dr_graph, DUMMY)
#     dr_graph.add_edges_from(connect_edges)
#
# # two_gate_profit = create_default_weight_func(deg_table_to_dict(deg_table))
# # profit_dict = dict(
# #     (tuple(sorted(npair)), two_gate_profit(*sorted(npair)))
# #     for npair in combinations(set(d_graph.nodes) | set([src, target]), 2)
# # )
#
#
# # @lru_cache(maxsize=67108864)
# # def profit_func(gates: tp.Iterable[Gate]) -> int:
# #     sorted_gates = [src] + sorted(gates) + [target]
# #     # ends_profits = 0
# #     # try:
# #     #     # Each gate is counted twice so we add 0.5 to account of each gate once
# #     #     ends_profits += sum(
# #     #         [
# #     #             profit_dict[(src, sorted_gates[0])],
# #     #             profit_dict[(sorted_gates[-1], target)],
# #     #         ]
# #     #     )
# #     # except IndexError:
# #     #     return np.inf
# #     profits = (profit_dict[gs] for gs in zip(sorted_gates[:-1], sorted_gates[1:]))
# #     return sum(profits)  # + len(gates)
#
# profit_dict = dict(
#     (gate, len(dna) - min(np.abs(gate.idx - var_pos_arr))) for gate in d_graph.nodes
# )
#
# for gate in d_graph.nodes:
#     if np.sign(gate.idx - var_pos_arr)[np.argmin(np.abs(gate.idx - var_pos_arr))] > 0:
#         profit_dict[gate] += 1
#     else:
#         profit_dict[gate] += 4
#
#
# @lru_cache(maxsize=67108864)
# def profit_func(gates: tp.Iterable[Gate]) -> int:
#     return sum(profit_dict[gate] for gate in gates)
#
#
# def gate_occupancy(gates: tp.Iterable[Gate]) -> tp.Tuple:
#     occ = [0] * (var_segemets.shape[0] + 1)
#     for gate in gates:
#         occ[var_segemets.shape[0] - sum(np.all(gate.idx < var_segemets, axis=1))] += 1
#     return tuple(occ)
#
#
# def valid_gates_func(var_pos_arr: np.array) -> tp.Callable[[tp.Iterable[Gate]], bool]:
#     # def valid_edges(gates: tp.Iterable[Gate]) -> bool:
#     #     if len(gates) < 2:
#     #         return True
#     #     gates = sorted(gates)
#     #     for g1, g2 in zip(gates[:-1], gates[1:]):
#     #         across_var_seg = np.any(
#     #             np.logical_and(g1.idx < var_pos_arr, var_pos_arr < g2.idx)
#     #         )
#     #         segment_len = g2.idx - g1.idx
#     #         if across_var_seg:
#     #             if segment_len < 20 or segment_len > 40:
#     #                 return False
#     #         else:
#     #             if segment_len < 20:
#     #                 return False
#     #     return True
#     # return valid_edges
#     def valid_edges(gates: tp.Iterable[Gate]) -> bool:
#         if len(gates) < 2:
#             return True
#         gates = sorted(gates)
#         avs = []
#         for g1, g2 in zip(gates[:-1], gates[1:]):
#             across_var_seg = np.any(
#                 np.logical_and(g1.idx < var_pos_arr, var_pos_arr < g2.idx)
#             )
#             segment_len = g2.idx - g1.idx
#             if across_var_seg:
#                 if segment_len < 20 or segment_len > 40:
#                     return False
#             else:
#                 if segment_len < 20:
#                     return False
#             avs.append(across_var_seg)
#         return any(avs)
#
#     return valid_edges
#
#
# max_occ = [2] * (var_segemets.shape[0] + 1)
# max_occ[0] = 1
# max_occ[-1] = 1
# max_occ = np.array(max_occ)
#
#
# def is_comp_occ(occ1: tp.Tuple, occ2: tp.Tuple) -> bool:
#     return not any(np.array(occ1) + np.array(occ2) > max_occ)
#
#
# edge_filter_func = valid_gates_func(var_pos_arr)
#
#
# def build_clique_tree(G: nx.Graph):
#     tree_struct, clique_dict, alpha = MCSCliqueTree(G)
#     rtree = tree_struct_to_tree(tree_struct, clique_dict)
#     root = [n for n, d in rtree.in_degree() if d == 0][0]
#
#     return rtree, root
#
#
# # for gs in combinations(dr_graph.nodes, 2):
# #     if not edge_filter_func(gs):
# #         dr_graph.add_edge(*gs)
#
# # H, alpha = nx.complete_to_chordal_graph(dr_graph)
#
# # profit_func = make_profit_func(
# #     deg_table_to_dict(deg_table),
# #     source=Gate(-1, "SOURCE"),
# #     target=Gate(len(dna), "TARGET"),
# # )
# # root = [n for n, d in rtree.in_degree() if d == 0][0]
# rtree, root = build_clique_tree(dr_graph)
# nx.drawing.nx_pydot.write_dot(rtree, "/home/labs/fleishman/shayho/Downloads/rtree.dot")
#
# f = nested_defaultdict(set, depth=3)
# f = algch8(rtree, root, f, profit_func, edge_filter_func, frozenset([DUMMY]))
# # f = algch9(rtree, root, f, profit_func)
#
# print(
#     sorted(
#         [ss for v_set in f[root] for d in f[root][v_set] for ss in f[root][v_set][d]],
#         key=lambda x: x[1],
#     )[:10]
# )
#
# var_edges = [edge for edge in d_graph.edges if is_var_edge(edge, var_pos_arr)]
# var_edges_props = dict(
#     (edge, d_graph.edges[edge])
#     for edge in d_graph.edges
#     if is_var_edge(edge, var_pos_arr)
# )
#
# var_edges_segments = [[] for i in range(var_segemets.shape[0])]
# for edge in var_edges:
#     i = np.argwhere(
#         np.any(
#             np.logical_and(edge[0].idx < var_segemets, var_segemets < edge[1].idx),
#             axis=1,
#         )
#     ).item()
#     var_edges_segments[i].append(edge)
#
# for i, var_seg in enumerate(var_edges_segments):
#     var_edges_segments[i] = sorted(var_seg, key=lambda x: var_edges_props[x]["weight"])
#
# take_top = 20
# sub_var_edges = list(
#     chain.from_iterable([var_seg[:take_top] for var_seg in var_edges_segments])
# )
#
# min_P = 0
# max_P = sum(
#     var_edges_props[x]["weight"]
#     for x in [var_seg[take_top - 1] for var_seg in var_edges_segments]
# )
#
# var_seg_g = nx.DiGraph()
# var_seg_g_edges = (
#     (vs1, vs2)
#     for ves1, ves2 in zip(var_edges_segments[:-1], var_edges_segments[1:])
#     for vs1, vs2 in product(ves1[:take_top], ves2[:take_top])
# )
# var_seg_g.add_edges_from(var_seg_g_edges)
#
# interval_g = nx.Graph()
# for edge1, edge2 in combinations(sub_var_edges, 2):
#     if edge1[1].idx < edge2[0].idx or edge2[1].idx < edge1[0].idx:
#         continue
#     interval_g.add_edge(edge1, edge2)
#
# bps_to_edge = defaultdict(set)
# for edge in sub_var_edges:
#     for u in edge:
#         bps_to_edge[u.bps].add(edge)
#
# sub_R_graph = R_graph.subgraph(bps_to_edge.keys()).copy()
# if not nx.is_connected(sub_R_graph):
#     DUMMY_BPS = "DUMMY"
#     connect_edges = create_connecet_edges(sub_R_graph, DUMMY_BPS)
#     sub_R_graph.add_edges_from(connect_edges)
# bps_to_edge[DUMMY_BPS] = frozenset()
#
# tw, Rtree = treewidth.treewidth_min_degree(sub_R_graph)
# tw2, Rtree2 = treewidth.treewidth_min_fill_in(sub_R_graph)
# if tw2 < tw:
#     Rtree = Rtree2
#
# bag_r_tree, node_bag_dict, max_id = baggify_tree(Rtree)
# root = find_suitable_root(bag_r_tree, DUMMY_BPS)
#
# bag_r_tree_bak = bag_r_tree.copy()
# root = make_nice(root, bag_r_tree, max_id + 1)
# rooted_bag_r_tree = nx.dfs_tree(bag_r_tree, root)
# for node in list(rooted_bag_r_tree.nodes):
#     if len(node.vertices) == 0:
#         rooted_bag_r_tree.remove_node(node)
# for node, in_deg in rooted_bag_r_tree.in_degree():
#     if in_deg == 0:
#         break
# old_root = root
# root = node
#
# nx.drawing.nx_pydot.write_dot(
#     rooted_bag_r_tree, "/home/labs/fleishman/shayho/Downloads/Rtree.dot"
# )
#
#
# bag_type = classify_bags(rooted_bag_r_tree, root)
# C = 12
#
#
# def profit_func(edge: tp.Tuple[Gate, Gate]):
#     return var_edges_props[edge]["weight"]
#
#
# def weight_func(edge: tp.Tuple[Gate, Gate]):
#     return C + 1 if DUMMY_BPS in segments_to_bps([edge]) else 1
#
#
# # try:
# #     with open("/home/labs/fleishman/shayho/Downloads/921_algTDCIG.pickle", "rb") as inf:
# #         sub_R_graph, rooted_bag_r_tree, root, interval_g, bps_to_edge, bag_type, profit_func, weight_func, min_P, max_P, C = pickle.load(
# #             inf
# #         )
# # except FileNotFoundError:
# #     with open(
# #         "/home/labs/fleishman/shayho/Downloads/921_algTDCIG.pickle", "wb"
# #     ) as outf:
# #         pickle.dump(
# #             [
# #                 sub_R_graph,
# #                 rooted_bag_r_tree,
# #                 root,
# #                 interval_g,
# #                 bps_to_edge,
# #                 bag_type,
# #                 profit_func,
# #                 weight_func,
# #                 min_P,
# #                 max_P,
# #                 C,
# #             ],
# #             outf,
# #         )
#
#
# f = AlgTDCIG(
#     sub_R_graph,
#     rooted_bag_r_tree,
#     root,
#     interval_g,
#     bps_to_edge,
#     bag_type,
#     profit_func,
#     weight_func,
#     min_P,
#     max_P,
#     C,
# )
# print(len(f))
#
# # min([(d, t[1]) for T in f[root][frozenset()].values() for d, t in T.items() if len(t[1]) == 3], key=lambda x: x[0])
#
# # TODO:
# # turn bag_r_tree to nice tree
# # root bag_r_tree to rooted tree
#
#
# # dna_seq = "CAAGTGCAATTACAAGAATCAGGGGGCGGGCTGGTCCAAGCCGGAGGTTCTCTGCGTTTATCCTGTGCGGCTTCTGGTTTCGACTTCTCTAAATACTACATGTCTTGGGTTCGTCAAGCCCCTGGTAAAGGCCTGGAATTCGTTGCGGCGATCTGGCCGTCTGGTTCTCACACCTACTACGCGGACTCTGTTAAAGGTCGCTTCACCATCTCTCGTGACAACTCTAAAAACCTGCTGTACCTGCAGATGAACTCTCTGCGTGCGGAAGACACCGCTGTATACTACTGCGCGCGTCACTCTCCGTACACCCACCAGCCGGACTACTGGGGACAGGGAACACAAGTAACTGTAAGTTCT"
# # positions = set(range(len(dna_seq)))
# # var_AA = {30, 31, 33, 56, 57, 59, 99, 101, 102, 104}
#
# # var_pos_list = [[3 * i, 3 * i + 1, 3 * i + 2] for i in var_AA]
# # var_pos = set(reduce(lambda x, y: x + y, var_pos_list))
# # var_pos_arr = np.array(sorted(var_pos))
# # const_pos = positions - var_pos
#
# # gm = GraphMaker(ggdata)
# # dna_g = make_default_graph(dna_seq, list(const_pos), list(var_pos), 20, 80, 20)
# # # for edge in G.edges():
# # #     G[edge[0]][edge[1]]["weight"] = edge_weight(edge, sub_var_pos_cost)
# # source_sink = [Gate(index=-1, bps="src"), Gate(index=len(dna_seq), bps="snk")]
# # incomp_graph.add_nodes_from(["src", "snk"])
#
# # bps_to_gate = defaultdict(set)
# # for edge in dna_g.edges:
# #     if is_var_edge(edge, var_pos_arr):
# #         bps_to_gate[edge[0].bps].add(edge[0])
# #         bps_to_gate[edge[1].bps].add(edge[1])
#
# # bps_count = dict((bps, len(gates)) for bps, gates in bps_to_gate.items())
# # bps_count["src"] = 1
# # bps_count["snk"] = 1
#
# # for edge in dna_g.edges:
# #     weight = 0.0
# #     if is_var_edge(edge, var_pos_arr):
# #         weight = count_edge_weight(edge, bps_count, incomp_graph)
# #     dna_g.edges[edge[0], edge[1]]["weight"] = weight
#
# # res = treewidth.treewidth_min_fill_in(G)
# # print(f"treewidth is {res[0]}")
# # tree = res[1]
#
# # # edges = [
# # #     (frozenset({"AACC", "GGTT"}), frozenset({"TGTT"})),
# # #     (frozenset({"AACC", "GGTT"}), frozenset({"ATTT"})),
# # #     (frozenset({"AACC", "GGTT"}), frozenset({"CTTT"})),
# # #     (frozenset({"AACC", "GGTT"}), frozenset({"GTTT"})),
# # #     (frozenset({"TGTT"}), frozenset({"AACA", "TGTT"})),
# # #     (frozenset({"ATTT"}), frozenset({"AAAT", "ATTT"})),
# # #     (frozenset({"CTTT"}), frozenset({"AAAG", "CTTT"})),
# # #     (frozenset({"GTTT"}), frozenset({"AAAC", "GTTT"})),
# # # ]
# # # tree = nx.Graph()
# # # tree.add_edges_from(edges)
# # nodes_id = dict((node, i) for i, node in enumerate(tree.node))
# # edges_ids = [(Bag(u, nodes_id[u]), Bag(v, nodes_id[v])) for u, v in tree.edges]
# # tree = nx.Graph()
# # tree.add_edges_from(edges_ids)
# # suitable_root = find_suitable_root(tree)
# # next_id = max(nodes_id.values()) + 1
# # nice_tree = tree.copy()
# # root = make_nice(suitable_root, nice_tree, next_id)
# # # # for i, node in enumerate(nice_tree.node):
# # # #     node.id = i
# # nx.drawing.nx_pydot.write_dot(
# #     nice_tree, "/home/labs/fleishman/shayho/Downloads/nice_tree.dot"
# # )
#
#
# # dna_seq = "CAAGTGCAATTACAAGAATCAGGGGGCGGGCTGGTCCAAGCCGGAGGTTCTCTGCGTTTATCCTGTGCGGCTTCTGGTTTCGACTTCTCTAAATACTACATGTCTTGGGTTCGTCAAGCCCCTGGTAAAGGCCTGGAATTCGTTGCGGCGATCTGGCCGTCTGGTTCTCACACCTACTACGCGGACTCTGTTAAAGGTCGCTTCACCATCTCTCGTGACAACTCTAAAAACCTGCTGTACCTGCAGATGAACTCTCTGCGTGCGGAAGACACCGCTGTATACTACTGCGCGCGTCACTCTCCGTACACCCACCAGCCGGACTACTGGGGACAGGGAACACAAGTAACTGTAAGTTCT"
# # positions = set(range(len(dna_seq)))
# # var_AA = {30, 31, 33, 56, 57, 59, 99, 101, 102, 104}
#
# # var_pos_list = [[3 * i, 3 * i + 1, 3 * i + 2] for i in var_AA]
# # var_pos = set(reduce(lambda x, y: x + y, var_pos_list))
# # const_pos = positions - var_pos
#
# # sub_dna = dna_seq[60:210]
# # sub_pos = set(range(len(dna_seq)))
# # sub_var_pos = set([i - 60 for i in var_pos if 60 <= i and i <= 210])
# # sub_var_pos_arr = np.array(sorted(sub_var_pos))
# # sub_const_pos = sub_pos - sub_var_pos
# # sub_var_pos_cost = {
# #     (30, 32): np.log(6),
# #     (33, 35): np.log(2),
# #     (39, 41): np.log(8),
# #     (108, 110): np.log(2),
# #     (111, 113): np.log(4),
# #     (117, 119): np.log(6),
# # }
#
#
# # ggdata = GGData()
# # gm = GraphMaker(ggdata)
# # G = gm.make_graph(dna_seq, list(const_pos), list(var_pos), 20, 80, 20)
# # # for edge in G.edges():
# # #     G[edge[0]][edge[1]]["weight"] = edge_weight(edge, sub_var_pos_cost)
# # source_sink = [Gate(index=-1, bps="src"), Gate(index=len(dna_seq), bps="snk")]
#
#
# # acc = defaultdict(list)
# # for edge in G.edges():
# #     cost = G.edges[edge[0], edge[1]]["weight"]
# #     if np.isinf(cost) or cost == 0:
# #         continue
# #     var_start_segmets = set()
# #     var_end_segments = set()
# #     for seg, _ in sub_var_pos_cost.items():
# #         if edge[0].index < seg[0]:
# #             var_start_segmets.add(seg[0])
# #         if edge[1].index > seg[1]:
# #             var_end_segments.add(seg[1])
# #     try:
# #         acc[(min(var_start_segmets), max(var_end_segments))].append((edge, cost))
# #     except ValueError:
# #         pass
#
#
# # restriction_g = nx.Graph()
# # restriction_g.add_nodes_from(ggdata.filter_self_binding_gates())
# # restriction_g.remove_nodes_from(source_sink)
#
# # edges = [
# #     (v1, v2)
# #     for v1, v2 in combinations(restriction_g.nodes, 2)
# #     if ggdata.gates_all_scores(v1.bps, v2.bps) >= 1000
# # ]
#
# # restriction_g.add_edges_from(edges)
# # connect_edges = create_connecet_edges(restriction_g, Gate(-1, "DUMMY"))
# # restriction_g.add_edges_from(connect_edges)
# # H, alpha = complete_to_chordal_graph(restriction_g)
# # alpha = LexBFSOrdering(restriction_g)
# # elim_order = [0] * len(alpha)
# # for key, val in alpha.items():
# #     elim_order[val - 1] = key
#
# # intersection_tree = intersection_graph(restriction_g, elim_order, alpha)
