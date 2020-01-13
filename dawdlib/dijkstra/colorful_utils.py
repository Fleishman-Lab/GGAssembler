import typing as tp
from collections import defaultdict
from functools import reduce

import networkx as nx
import numpy as np
from dawdlib.dijkstra import colorful
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gate_restriction import gen_gate_restriction_graph
from networkx.algorithms.shortest_paths.weighted import _weight_function


def _MCSCliqueTree(
    G: nx.Graph
) -> tp.Tuple[tp.List[tp.Tuple[int, int]], tp.Dict[int, tp.Set], tp.Dict[tp.Any, int]]:
    """
    Ref: Galinier et al. - Chordal graphs and their clique graphs (2005)
    Data: A graph G = (V, E)
    Result: If the input graph is chordal: a PEO and an associated
    clique-tree T = (I, F) where I is the set of maximal cliques
    begin
        each vertex of X is initialized with the empty set
        previousmark = -1
        j=0
        for i=n to 1 do
            choose a vertex x not yet numbered such that |mark(x)| is maximum
            if mark(x) < previousmark then
                j=j+1
                create the maximal clique Cj = M(x) U {x}
                create the tie between Cj and C(last(x))
            else
                Cj = Cj U {x}
            for each y neighbour of x do
                M(y) = M(y) U {x}
                mark(y) = mark(y) + 1
                last(y) = x
            previousmark = mark(x)
            x is numbered by i
            C(x) = j
    end
    """
    heap = dict((n, 0) for n in G.nodes)
    T: tp.List[tp.Tuple[int, int]] = []
    alpha: tp.Dict[tp.Any, int] = dict()
    C: tp.Dict[int, tp.Set] = defaultdict(set)
    C_dict: tp.Dict[tp.Any, int] = dict((n, 0) for n in G.nodes)
    M: tp.Dict[tp.Any, tp.Set] = dict((n, set()) for n in G.nodes)
    previousmark = -1
    j = 0
    last: tp.Dict[tp.Any, int] = dict((n, 0) for n in G.nodes)
    for i in range(len(G.nodes), 0, -1):
        u = max(heap, key=heap.get)
        marku = heap.pop(u)
        if marku <= previousmark:
            j += 1
            C[j] = M[u] | set([u])
            T.append((j, C_dict[last[u]]))
        else:
            C[j] = C[j] | set([u])
        for v in G[u]:
            try:
                heap[v] += 1
                M[v] = M[v] | set([u])
                last[v] = u
            except KeyError:
                pass
        previousmark = marku
        alpha[u] = i
        C_dict[u] = j
    return T, C, alpha


def _enumerate_cliques(graph: nx.Graph) -> tp.Dict[int, tp.Set[tp.Any]]:
    cliques: tp.List[tp.Set] = []
    for cc in nx.connected_components(graph):
        tree_struct, clique_dict, alpha = _MCSCliqueTree(graph.subgraph(cc))
        cliques.extend(clique_dict.values())

    return dict(enumerate(cliques))


def _color_vertices(
    clique_colors: tp.Dict[int, tp.Set[tp.Any]]
) -> tp.DefaultDict[tp.Any, tp.FrozenSet[int]]:
    vertex_colors = defaultdict(frozenset)
    for color, clique in clique_colors.items():
        for vertex in clique:
            vertex_colors[vertex] |= frozenset([color])
    return vertex_colors


def _color_gates(
    graph: tp.Union[nx.Graph, nx.DiGraph],
    vertex_colors: tp.DefaultDict[tp.Any, tp.FrozenSet[int]],
) -> tp.Dict[tp.Any, tp.FrozenSet[int]]:
    gate_colors = {}
    for node in graph.nodes:
        gate_colors[node] = vertex_colors[node.bps]
    return gate_colors


def color_gates(
    graph: tp.Union[nx.Graph, nx.DiGraph],
    ggdata: GGData,
    gate_self_binding_min: int = 2000,
    gate_crosstalk_max: int = 1000,
) -> tp.Dict[tp.Any, tp.FrozenSet[int]]:
    r_graph = gen_gate_restriction_graph(
        ggdata, gate_self_binding_min, gate_crosstalk_max
    )
    clique_colors = _enumerate_cliques(r_graph)
    vertex_colors = _color_vertices(clique_colors)
    return _color_gates(graph, vertex_colors)


def _nxweighttonp(
    graph: tp.Union[nx.Graph, nx.DiGraph], nodes: tp.List, weight="weight"
) -> np.array:
    weight_fn = _weight_function(graph, weight)
    weight_arr = np.empty((len(nodes), len(nodes)), dtype=np.uint16)
    weight_arr[:] = np.iinfo(np.uint16).max
    for edge in graph.edges:
        weight_arr[nodes.index(edge[0]), nodes.index(edge[1])] = weight_fn(
            edge[0], edge[1], graph.edges[edge[0], edge[1]]
        )
    return weight_arr


def _nxcolortonp(color: tp.Callable[[tp.Any], tp.FrozenSet[int]], nodes: tp.List):
    colors = np.unique(list(reduce(set.union, [set(color(node)) for node in nodes])))
    color_arr = np.zeros((len(nodes), colors.shape[0]), dtype=np.ubyte)
    for node in nodes:
        n_index = nodes.index(node)
        for node_color in color(node):
            color_arr[n_index, np.flatnonzero(colors == node_color).item()] = True
    return color_arr


def _nxgraphtomap(weight_arr: np.array) -> tp.Dict:
    cgraph = defaultdict(list)
    rows, cols = np.nonzero(weight_arr < np.iinfo(np.uint16).max)
    for i, j in zip(rows, cols):
        cgraph[int(i)].append(int(j))
    cgraph = dict(cgraph)

    return colorful.graphdict2map(cgraph)


def nxtonumpy(graph, sources, target, color, no_colors, weight="weight"):
    nodes = sorted(graph.nodes)

    weight_arr = _nxweighttonp(graph, nodes, weight)
    color_arr = _nxcolortonp(color, nodes)
    cgraph = _nxgraphtomap(weight_arr)
    sources = np.array([nodes.index(source) for source in sources], dtype=np.int32)
    target = nodes.index(target)
    new_color = np.zeros((color_arr.shape[0],), dtype=np.uint16)
    pred = np.zeros_like(weight_arr, dtype=np.ubyte)
    dist = np.empty((weight_arr.shape[0], 2 ** no_colors), np.uint16)
    seen = np.empty_like(dist)

    return (
        nodes,
        cgraph,
        weight_arr,
        sources,
        target,
        color_arr,
        new_color,
        no_colors,
        pred,
        dist,
        seen,
    )
