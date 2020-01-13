# import typing as tp
# import networkx as nx
# import numpy as np
#
# from itertools import combinations, product
# from functools import reduce
# from dawdlib.golden_gate.gate_data import GGData
# from dawdlib.golden_gate.graph_maker import GraphMaker
# from dawdlib.golden_gate.gate import Gate
#
#
# from pymoo.model.problem import Problem
#
# class MyProblem(Problem):
#
#     def __init__(self, n_var, n_obj = )
#
# def find_cliques(G, source_n: set, target_n: set, min_depth: int):
#     """Returns all maximal cliques in an undirected graph.
#
#     For each node *v*, a *maximal clique for v* is a largest complete
#     subgraph containing *v*. The largest maximal clique is sometimes
#     called the *maximum clique*.
#
#     This function returns an iterator over cliques, each of which is a
#     list of nodes. It is an iterative implementation, so should not
#     suffer from recursion depth issues.
#
#     Parameters
#     ----------
#     G : NetworkX graph
#         An undirected graph.
#
#     Returns
#     -------
#     iterator
#         An iterator over maximal cliques, each of which is a list of
#         nodes in `G`. The order of cliques is arbitrary.
#
#     See Also
#     --------
#     find_cliques_recursive
#         A recursive version of the same algorithm.
#
#     Notes
#     -----
#     To obtain a list of all maximal cliques, use
#     `list(find_cliques(G))`. However, be aware that in the worst-case,
#     the length of this list can be exponential in the number of nodes in
#     the graph. This function avoids storing all cliques in memory by
#     only keeping current candidate node lists in memory during its search.
#
#     This implementation is based on the algorithm published by Bron and
#     Kerbosch (1973) [1]_, as adapted by Tomita, Tanaka and Takahashi
#     (2006) [2]_ and discussed in Cazals and Karande (2008) [3]_. It
#     essentially unrolls the recursion used in the references to avoid
#     issues of recursion stack depth (for a recursive implementation, see
#     :func:`find_cliques_recursive`).
#
#     This algorithm ignores self-loops and parallel edges, since cliques
#     are not conventionally defined with such edges.
#
#     References
#     ----------
#     .. [1] Bron, C. and Kerbosch, J.
#        "Algorithm 457: finding all cliques of an undirected graph".
#        *Communications of the ACM* 16, 9 (Sep. 1973), 575--577.
#        <http://portal.acm.org/citation.cfm?doid=362342.362367>
#
#     .. [2] Etsuji Tomita, Akira Tanaka, Haruhisa Takahashi,
#        "The worst-case time complexity for generating all maximal
#        cliques and computational experiments",
#        *Theoretical Computer Science*, Volume 363, Issue 1,
#        Computing and Combinatorics,
#        10th Annual International Conference on
#        Computing and Combinatorics (COCOON 2004), 25 October 2006, Pages 28--42
#        <https://doi.org/10.1016/j.tcs.2006.06.015>
#
#     .. [3] F. Cazals, C. Karande,
#        "A note on the problem of reporting maximal cliques",
#        *Theoretical Computer Science*,
#        Volume 407, Issues 1--3, 6 November 2008, Pages 564--568,
#        <https://doi.org/10.1016/j.tcs.2008.05.010>
#
#     """
#     if len(G) == 0:
#         return
#
#     adj = {u: {v for v in G[u] if v != u} for u in G}
#     Q = [None]
#
#     subg = set(G)
#     cand = set(G)
#     u = max(subg, key=lambda u: len(cand & adj[u]))
#     ext_u = cand - adj[u]
#     stack = []
#
#     try:
#         while True:
#             if ext_u:
#                 q = ext_u.pop()
#                 cand.remove(q)
#                 Q[-1] = q
#                 adj_q = adj[q]
#                 subg_q = subg & adj_q
#                 if (
#                     not subg_q
#                     and len(Q) > min_depth
#                     and len(source_n.intersection(Q)) > 0
#                     and len(target_n.intersection(Q)) > 0
#                 ):
#                 # if not subg_q and len(Q) >= min_depth:
#                     yield Q[:]
#                 else:
#                     cand_q = cand & adj_q
#                     if cand_q:
#                         stack.append((subg, cand, ext_u))
#                         Q.append(None)
#                         subg = subg_q
#                         cand = cand_q
#                         u = max(subg, key=lambda u: len(cand & adj[u]))
#                         ext_u = cand - adj[u]
#             else:
#                 Q.pop()
#                 subg, cand, ext_u = stack.pop()
#     except IndexError:
#         pass
#
#
# def find_cliques_recursive(
#     G: nx.Graph, source_n: set, target_n: set, min_depth: int, max_depth: int
# ) -> tp.Generator[tp.List, None, None]:
#     """Reimplementation of find_cliques_recursive of networkx.
#
#     This implementation requires that the cliques includes the source and target and are of size at least min_depth.
#     It starts from the vertex with the least amount of neighbours.
#
#     The recursion is limited in depth by max depth.
#
#     Arguments:
#         G {[type]} -- [description]
#         source {[type]} -- [description]
#         target {[type]} -- [description]
#         min_depth {[type]} -- [description]
#         max_depth {[type]} -- [description]
#     """
#
#     if len(G) == 0:
#         return iter([])
#
#     adj = {u: {v for v in G[u] if v != u} for u in G}
#     Q = []
#
#     def expand(subg, cand, depth):
#         u = max(subg, key=lambda u: len(cand & adj[u]))
#         for q in cand - adj[u]:
#             cand.remove(q)
#             Q.append(q)
#             adj_q = adj[q]
#             subg_q = subg & adj_q
#             if (
#                 not subg_q
#                 and depth > min_depth
#                 and len(source_n.intersection(Q))
#                 and len(target_n.intersection(Q))
#             ):
#                 yield Q[:]
#             else:
#                 cand_q = cand & adj_q
#                 if cand_q:
#                     for clique in expand(subg_q, cand_q, depth + 1):
#                         yield clique
#             Q.pop()
#
#     return expand(set(G), set(G), 0)
#
#
# def verify_clique(G: nx.Graph, nodes: set) -> bool:
#     for u, v in combinations(nodes, 2):
#         if not G.has_edge(u, v):
#             return False
#     return True
#     # clique = [set(G[u]) for u in nodes]
#     # clique = reduce(lambda u, v: u & v, clique)
#     # return nodes.issubset(clique)
#
#
# def k_cliques(graph: nx.Graph, source, target, min_depth, max_depth):
#     # 2-cliques
#     cliques = [{source, j} for j in graph[source] if source != j]
#     k = 2
#
#     while cliques:
#         # result
#         yield k, list(filter(lambda clique: target in clique, cliques))
#         k += 1
#         if k > max_depth:
#             break
#
#         print(f"Doing cliques of size {k}.")
#         # merge k-cliques into (k+1)-cliques
#         cliques_1 = set()
#         for u, v in combinations(cliques, 2):
#             w = u ^ v
#             if len(w) == 2 and graph.has_edge(*w):
#                 cliques_1.add(tuple(u | w))  # remove duplicates
#         cliques = list(map(set, cliques_1))
#
#
# class FakeGGData():
#
#     def __init__(self, gates):
#         self.gates = gates
#
#     def gates_all_scores(self, *args):
#         return 0
#
#     def filter_self_binding_gates(self, *args):
#         return self.gates
#
#
#
# ggdata = GGData()
# nodes = ggdata.filter_self_binding_gates()
# # src_sink = ["src", "snk"]
# edges = [
#     (v1, v2)
#     for v1, v2 in combinations(nodes, 2)
#     if ggdata.gates_all_scores(v1, v2) >= 1000
# ]
# # edges += [(v1, v2) for v1, v2 in product(src_sink, nodes)]
# # edges += [("src", "snk")]
# # nodes += src_sink
# # G = nx.Graph()
# # G.add_nodes_from(nodes)
# # G.add_edges_from(edges)
# dna_seq = "CAAGTGCAATTACAAGAATCAGGGGGCGGGCTGGTCCAAGCCGGAGGTTCTCTGCGTTTATCCTGTGCGGCTTCTGGTTTCGACTTCTCTAAATACTACATGTCTTGGGTTCGTCAAGCCCCTGGTAAAGGCCTGGAATTCGTTGCGGCGATCTGGCCGTCTGGTTCTCACACCTACTACGCGGACTCTGTTAAAGGTCGCTTCACCATCTCTCGTGACAACTCTAAAAACCTGCTGTACCTGCAGATGAACTCTCTGCGTGCGGAAGACACCGCTGTATACTACTGCGCGCGTCACTCTCCGTACACCCACCAGCCGGACTACTGGGGACAGGGAACACAAGTAACTGTAAGTTCT"
# positions = set(range(len(dna_seq)))
# var_AA = {30, 31, 33, 56, 57, 59, 99, 101, 102, 104}
# var_pos = [[3 * i, 3 * i + 1, 3 * i + 2] for i in var_AA]
# var_pos = set(reduce(lambda x, y: x + y, var_pos))
# const_pos = positions - var_pos
# ggdata = FakeGGData(['ACAG'])
# gm = GraphMaker(ggdata)
# G = gm.make_graph(dna_seq, list(const_pos), list(var_pos), 20, 80, 20)
# G = G.to_undirected()
# # src_sink = [Gate(index=-1, bps="src"), Gate(index=357, bps="snk")]
# # edges = [(v1, v2) for v1, v2 in product(src_sink, G.nodes)]
# # G.add_edges_from(edges)
# # source_n = set(G[Gate(index=-1, bps="src")])
# # target_n = set(G[Gate(index=357, bps="snk")])
# # G.remove_nodes_from(src_sink)
# # np.random.seed(0)
# # mat = np.zeros((300, 300))
# # ones = np.random.choice(300 ** 2, int(np.ceil(0.14 * 300 ** 2)))
# # mat.reshape(-1)[ones] = 1
# # G = nx.Graph(mat)
#
# # source_n = set(G[0])
# # target_n = set(G[299])
# # G.remove_nodes_from([0, 299])
#
# # count = 0
# # # for sol in find_cliques(G, source_n, target_n, 4):
# # #     # print(verify_clique(G, set(sol)), sol)
# # #     count += 1
#
# # for sol in nx.find_cliques(G):
# #     count += 1
#
# # print(count)
#
#
# # print(len(res), res)
# # for sol in find_cliques_recursive(G, "src", "snk", -1, 10):
# #     print(verify_clique(G, set(sol)), sol)
#
# # for sol in k_cliques(G, "src", "snk", -1, 10):
# #     print(sol)
