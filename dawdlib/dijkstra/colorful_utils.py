import typing as tp
from collections import OrderedDict, defaultdict
from functools import reduce

import networkx as nx
import numpy as np
from networkx.algorithms.shortest_paths.weighted import _weight_function

from dawdlib.dijkstra import colorful
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gate_restriction import gen_gate_restriction_graph


def _MCSCliqueTree(
    G: nx.Graph,
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
            C[j] = M[u] | {u}
            T.append((j, C_dict[last[u]]))
        else:
            C[j] |= {u}
        for v in G[u]:
            try:
                heap[v] += 1
                M[v] |= {u}
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
        _, clique_dict, _ = _MCSCliqueTree(graph.subgraph(cc))
        cliques.extend(clique_dict.values())

    return dict((i, c) for i, c in enumerate(cliques))


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
    r_graph: nx.Graph = None,
) -> tp.Tuple[tp.Dict[tp.Any, tp.FrozenSet[int]], tp.List[int]]:
    if r_graph is None:
        r_graph = gen_gate_restriction_graph(
            ggdata, gate_self_binding_min, gate_crosstalk_max
        )
    clique_colors = _enumerate_cliques(r_graph)
    vertex_colors = _color_vertices(clique_colors)
    return _color_gates(graph, vertex_colors), num_of_colors(r_graph, vertex_colors)


def num_of_colors(
    r_graph: nx.Graph, vertex_colors: tp.Dict[tp.Any, tp.FrozenSet[int]]
) -> tp.List[int]:
    num_colors: tp.List[int] = []
    vertices: tp.List[tp.Any] = []
    vertex_colors = OrderedDict(
        sorted(vertex_colors.items(), key=lambda x: len(x[1]), reverse=True)
    )
    for key, val in vertex_colors.items():
        if not r_graph.subgraph(vertices + [key]).number_of_edges():
            vertices.append(key)
            num_colors.append(len(val))
    return num_colors


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


def _nxcolortonp(
    color: tp.Callable[[tp.Any], tp.FrozenSet[int]], nodes: tp.List
) -> np.array:
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


def nxtonumpy(
    graph: nx.Graph,
    sources: tp.List,
    target,
    color: tp.Callable[[tp.Any], tp.FrozenSet[int]],
    no_colors: int,
    weight="weight",
):
    nodes = list(graph.nodes)

    l_nodes = np.log2(len(nodes))
    l_nodes = l_nodes if l_nodes > no_colors else no_colors

    weight_arr = _nxweighttonp(graph, nodes, weight)
    color_arr = _nxcolortonp(color, nodes)
    cgraph = _nxgraphtomap(weight_arr)
    sources = np.array([nodes.index(source) for source in sources], dtype=np.int32)
    target = nodes.index(target)
    pred = np.zeros_like(weight_arr, dtype=np.ubyte)

    dtype = np.uint8
    if l_nodes < 9:
        dtype = np.uint8
    elif 8 < l_nodes < 17:
        dtype = np.uint16
    elif 16 < l_nodes < 33:
        dtype = np.uint32
    elif 33 < l_nodes < 65:
        dtype = np.uint64
    else:
        raise ValueError("The required number of colors exceeds 64. Please adjust!")

    new_color = np.zeros((color_arr.shape[0],), dtype=dtype)
    dist = np.empty((weight_arr.shape[0], 2 ** no_colors), dtype=dtype)
    seen = np.empty_like(dist)
    s_pred = np.empty_like(dist)

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
        s_pred,
        dist,
        seen,
    )


def yield_shortest_paths(nodes, src, trgt, clr_map, pred, limit=0):
    pred_dict = defaultdict(list)
    for key, val in np.argwhere(pred).tolist():
        pred_dict[key].append(val)
    stack = [[trgt, 0, clr_map[trgt], 0]]
    top = 0
    while top >= 0:
        node, i, p_clr, p_len = stack[top]
        if node == src:
            yield [nodes[p] for p, n, c, l in reversed(stack[: top + 1])]
        if len(pred_dict[node]) > i:
            p_node = pred_dict[node][i][0]
            if not (p_clr & clr_map[p_node]) and not 0 < limit < p_len:
                top += 1
                p_tup = [
                    p_node,
                    0,
                    p_clr | clr_map[p_node],
                    p_len + 1,
                ]
                if top == len(stack):
                    stack.append(p_tup)
                else:
                    stack[top] = p_tup
            else:
                stack[top][1] += 1
        else:
            stack[top - 1][1] += 1
            top -= 1


def return_shortest_path(nodes, src, trgt, s_pred, seen, clr_map, limit=0):
    stack = [[trgt, np.argmin(seen[trgt]), 0]]
    top = 0
    while top >= 0:
        v, p_col, p_len = stack[top]
        if v == src:
            return [nodes[p] for p, _, _ in reversed(stack[: top + 1])]
        if 0 < limit < p_len:
            raise nx.NetworkXNoPath()
        u = s_pred[v, p_col]
        top += 1
        u_tup = [u, p_col - clr_map[v], p_len + 1]
        stack.append(u_tup)
