import typing as tp
from collections import OrderedDict, defaultdict
from functools import reduce
from itertools import chain

import networkx as nx
import numpy as np
from networkx.algorithms.shortest_paths.weighted import _weight_function

# from dawdlib.dijkstra import colorful_algorithm  # pytype: disable=import-error
from dawdlib.dijkstra import colourful_dijkstra  # pytype: disable=import-error
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gate_restriction import gen_gate_restriction_graph


class ShortestPathFinder:
    def __init__(self, graph: nx.Graph, ggdata: GGData, source: tp.Any, target: tp.Any):
        self.graph = graph
        self.ggdata = ggdata
        self.source = source
        self.target = target

        self.gate_colors, self.num_of_colors = color_gates(self.graph, self.ggdata)
        self.node_map = dict((v, k) for k, v in enumerate(graph.nodes))
        self.idx_node_map = dict((k, v) for k, v in enumerate(graph.nodes))
        self.graph_edges = [
            (self.node_map[u], self.node_map[v], c)
            if v != source
            else (self.node_map[v], self.node_map[u], c)
            for (u, v, c) in graph.edges.data("weight", default=0)
        ]

        self.all_gate_colors = set()
        for val in self.gate_colors.values():
            self.all_gate_colors.update(val)

    def random_recolor(
        self, len_cutoff: tp.Optional[int], no_colors: tp.Optional[int]
    ) -> tp.Dict[int, int]:
        if no_colors is None and len_cutoff is None:
            raise ValueError(
                "At least one of len_cutoff or no_colors must be provided."
            )
        if no_colors is None:
            no_colors = sum(self.num_of_colors[:len_cutoff])
        dtype = np.uint8
        if no_colors > 8:
            dtype = np.uint16
        if no_colors > 16:
            dtype = np.uint32
        if no_colors > 32:
            dtype = np.uint64
        if no_colors > 64:
            raise ValueError(
                f"Len limit of {len_cutoff} results in over 64 colors please provide limit on the number of colors."
            )
        selected_colors = (2 ** np.arange(no_colors)).astype(dtype)
        color_mapping = dict(
            (k, np.random.choice(selected_colors).astype(dtype))
            for k in self.all_gate_colors
        )
        return dict(
            (
                self.node_map[n],
                np.sum(
                    [color_mapping[c] for c in self.gate_colors[n]], dtype=dtype
                ).item(),
            )
            for n in self.graph.nodes
        )

    def find_shortest_path(
        self, len_cutoff: tp.Optional[int] = None, no_colors: tp.Optional[int] = None
    ):
        if no_colors is None and len_cutoff is None:
            raise ValueError(
                "At least one of len_cutoff or no_colors must be provided."
            )
        path = colourful_dijkstra.colourful_shortest_path(
            self.graph_edges,
            self.node_map[self.source],
            self.node_map[self.target],
            self.random_recolor(len_cutoff, no_colors),
            len_cutoff,
        )
        return [self.idx_node_map[idx] for idx in path]


def all_shortest_paths(
    ggdata,
    G,
    source,
    target,
    no_colors=0,
    len_cutoff=0,
    weight=None,
    method="dijkstra",
    retries=1,
    res_graph: nx.Graph = None,
):
    """Compute all shortest paths in the graph.

    Parameters
    ----------
    ggdata: GGData
       GGdata for the reaction.

    G : NetworkX graph

    source : node
       Starting node for path.

    target : node
       Ending node for path.

    no_colors: int
        The number of colors to use in the random recoloring.
        The currently supported maximum is set to 18.

    len_cutoff: int, optional (default = 0)
        The maximum path length (numbe of gates in the solution)

    weight : None or string, optional (default = None)
       If None, every edge has weight/distance/cost 1.
       If a string, use this edge attribute as the edge weight.
       Any edge attribute not present defaults to 1.

    method : string, optional (default = 'dijkstra')
       The algorithm to use to compute the path lengths.
       Supported options: 'dijkstra'.
       Other inputs produce a ValueError.

    retries: int, optional (default = 1)
        The number of random coloring retries to check before raising NetworkXNoPath exception

    res_graph: nx.Graph, optional (default = None)
        Allows the user to provide a custom restriction graph used to color nodes

    Returns
    -------
    paths : generator of lists
        A generator of all paths between source and target.

    Raises
    ------
    ValueError
        If `method` is not among the supported options.
        If both `no_colors` and `len_cutoff` are not provided.

    NetworkXNoPath
        If `target` cannot be reached from `source`.


    Notes
    -----
    There may be many shortest paths between the source and target.

    """
    (
        nodes,
        cgraph,
        weight_arr,
        sources,
        trgt,
        color_arr,
        new_color,
        no_colors,
        pred,
        s_pred,
        dist,
        seen,
    ) = _prep_data(
        ggdata,
        G,
        source,
        target,
        no_colors=no_colors,
        len_cutoff=len_cutoff,
        weight=weight,
        method=method,
        retries=retries,
        res_graph=res_graph,
    )
    return all_colorful_shortest_paths(
        nodes,
        cgraph,
        weight_arr,
        sources,
        trgt,
        color_arr,
        new_color,
        no_colors,
        pred,
        s_pred,
        dist,
        seen,
        limit=len_cutoff,
    )


def shortest_path(
    ggdata,
    G,
    source,
    target,
    no_colors=0,
    len_cutoff=0,
    weight=None,
    method="dijkstra",
    retries=1,
    res_graph: nx.Graph = None,
):
    """Compute all shortest paths in the graph.

    Parameters
    ----------
    ggdata: GGData
       GGdata for the reaction.

    G : NetworkX graph

    source : node
       Starting node for path.

    target : node
       Ending node for path.

    no_colors: int
        The number of colors to use in the random recoloring.
        The currently supported maximum is set to 18.

    len_cutoff: int, optional (default = 0)
        The maximum path length

    weight : None or string, optional (default = None)
       If None, every edge has weight/distance/cost 1.
       If a string, use this edge attribute as the edge weight.
       Any edge attribute not present defaults to 1.

    method : string, optional (default = 'dijkstra')
       The algorithm to use to compute the path lengths.
       Supported options: 'dijkstra'.
       Other inputs produce a ValueError.

    retries: int, optional (default = 1)
        The number of random coloring retries to check before raising NetworkXNoPath exception

    res_graph: nx.Graph, optional (default = None)
        Allows the user to provide a custom restriction graph used to color nodes

    Returns
    -------
    paths : generator of lists
        A generator of all paths between source and target.

    Raises
    ------
    ValueError
        If `method` is not among the supported options.
        If both `no_colors` and `len_cutoff` are not provided.

    NetworkXNoPath
        If `target` cannot be reached from `source`.


    Notes
    -----
    There may be many shortest paths between the source and target.

    """
    (
        nodes,
        cgraph,
        weight_arr,
        sources,
        trgt,
        color_arr,
        new_color,
        no_colors,
        pred,
        s_pred,
        dist,
        seen,
    ) = _prep_data(
        ggdata,
        G,
        source,
        target,
        no_colors=no_colors,
        len_cutoff=len_cutoff,
        weight=weight,
        method=method,
        retries=retries,
        res_graph=res_graph,
    )
    return colorful_shortest_path(
        nodes,
        cgraph,
        weight_arr,
        sources,
        trgt,
        color_arr,
        new_color,
        no_colors,
        pred,
        s_pred,
        dist,
        seen,
        limit=len_cutoff,
    )


def _prep_data(
    ggdata,
    G,
    source,
    target,
    no_colors=0,
    len_cutoff=0,
    weight=None,
    method="dijkstra",
    retries=1,
    res_graph: nx.Graph = None,
):
    if method != "dijkstra":
        raise ValueError("method not supported: {}".format(method))
    if not no_colors and not len_cutoff:
        raise ValueError("Either no_colors or len_cutoff are required.")

    gate_colors, num_of_colors = color_gates(G, ggdata, r_graph=res_graph)
    if not no_colors:
        no_colors = sum(num_of_colors[:len_cutoff])
    if no_colors > 18:
        print(
            f"Warning! number of colors required {no_colors} is larger than recommended."
        )
    return nxtonumpy(G, [source], target, gate_colors.get, no_colors)


def all_colorful_shortest_paths(
    nodes: tp.List[int],
    G: tp.Dict[int, tp.FrozenSet[int]],
    weight: np.ndarray,
    sources: np.ndarray,
    target: int,
    colors: np.ndarray,
    color_map: np.ndarray,
    no_colors: int,
    pred: np.ndarray,
    s_pred: np.ndarray,
    dist: np.ndarray,
    seen: np.ndarray,
    limit: int = 0,
    repetitions: int = 1,
) -> tp.Iterable[tp.Any]:
    """
    :param list nodes: The list of nodes, used to identify each node with an index
    :param dict G: A mapping representation of the golden gate graph
    :param numpy.ndarray weight: A numpy array representations of graph edge weights
    :param numpy.ndarray sources: A numpy array holding indices of sources
    :param int target: The path target
    :param numpy.ndarray colors: A numpy array mapping between a node and it's colors
    :param numpy.ndarray color_map: A numpy array holding the random recoloring
    :param int no_colors: How many colors were used to created the dist and seen arrays
    :param numpy.ndarray pred: A numpy array pointing to the predecessor of each node in the paths found
    :param numpy.ndarray s_pred: A numpy array pointing to a single predecessor of each node
    :param numpy.ndarray dist: A numpy array holding the minimum distance from source to node using specific colors
    :param numpy.ndarray seen: A numpy array holding the minimum distance from source to node using specific colors
    :param int limit: A limit to the number of gates found in the shortest path
    :param int repetitions: How many times to try and find a shortest path using a random recoloring of no_colors
    :return: A generator yielding shortest colorful paths
    :rtype: generator
    """
    found, res = colorful_algorithm.find_shortest_paths(
        G,
        weight,
        sources,
        colors,
        color_map,
        no_colors,
        pred,
        s_pred,
        dist,
        seen,
        limit,
        repetitions,
    )
    if found:
        return chain.from_iterable(
            (
                yield_shortest_paths(nodes, src, target, color_map, pred)
                for src in sources
            )
        )
    else:
        raise nx.NetworkXNoPath


def colorful_shortest_path(
    nodes: tp.List[int],
    G: tp.Dict[int, tp.FrozenSet[int]],
    weight: np.ndarray,
    sources: np.ndarray,
    target: int,
    colors: np.ndarray,
    color_map: np.ndarray,
    no_colors: int,
    pred: np.ndarray,
    s_pred: np.ndarray,
    dist: np.ndarray,
    seen: np.ndarray,
    limit: int = 0,
    repetitions: int = 1,
) -> tp.Iterable[tp.Any]:
    """
    :param list nodes: The list of nodes, used to identify each node with an index
    :param dict G: A mapping representation of the golden gate graph
    :param numpy.ndarray weight: A numpy array representations of graph edge weights
    :param numpy.ndarray sources: A numpy array holding indices of sources
    :param int target: The path target
    :param numpy.ndarray colors: A numpy array mapping between a node and it's colors
    :param numpy.ndarray color_map: A numpy array holding the random recoloring
    :param int no_colors: How many colors were used to created the dist and seen arrays
    :param numpy.ndarray pred: A numpy array pointing to the predecessor of each node in the paths found
    :param numpy.ndarray s_pred: A numpy array pointing to a single predecessor of each node
    :param numpy.ndarray dist: A numpy array holding the minimum distance from source to node using specific colors
    :param numpy.ndarray seen: A numpy array holding the minimum distance from source to node using specific colors
    :param int limit: A limit to the number of gates found in the shortest path
    :param int repetitions: How many times to try and find a shortest path using a random recoloring of no_colors
    :return: A generator yielding shortest colorful paths
    :rtype: generator
    """
    found, res = colorful_algorithm.find_shortest_paths(
        G,
        weight,
        sources,
        colors,
        color_map,
        no_colors,
        pred,
        s_pred,
        dist,
        seen,
        limit,
        repetitions,
    )
    if found:
        return chain(
            return_shortest_path(nodes, src, target, s_pred, seen, color_map, limit)
            for src in sources
        )
    else:
        raise nx.NetworkXNoPath


def _color_gates(
    graph: tp.Union[nx.Graph, nx.DiGraph],
    vertex_colors: nx.Graph,
) -> tp.Dict[tp.Any, tp.FrozenSet[int]]:
    gate_colors = {}
    for node in graph.nodes:
        gate_colors[node] = (
            frozenset(vertex_colors[node.bps])
            if not node.src_or_target
            else frozenset()
        )
    return gate_colors


def color_gates(
    graph: tp.Union[nx.Graph, nx.DiGraph],
    ggdata: GGData,
    r_graph: nx.Graph = None,
) -> tp.Tuple[tp.Dict[tp.Any, tp.FrozenSet[int]], tp.List[int]]:
    if r_graph is None:
        r_graph = gen_gate_restriction_graph(ggdata)
    graph_bps = set([g.bps for g in graph.nodes])
    subr_graph = r_graph.subgraph(graph_bps)
    vertex_cliques = nx.make_clique_bipartite(subr_graph)
    return _color_gates(graph, vertex_cliques), num_of_colors(
        subr_graph, vertex_cliques
    )


def num_of_colors(r_graph: nx.Graph, vertex_cliques: nx.Graph) -> tp.List[int]:
    num_colors: tp.List[int] = []
    vertices: tp.List[tp.Any] = []
    vertex_colors = dict((node, len(vertex_cliques[node])) for node in r_graph.nodes)
    vertex_colors = OrderedDict(
        sorted(vertex_colors.items(), key=lambda x: x[1], reverse=True)
    )
    for key, val in vertex_colors.items():
        if not r_graph.subgraph(vertices + [key]).number_of_edges():
            vertices.append(key)
            num_colors.append(val)
    return num_colors


def _nxweighttonp(
    graph: tp.Union[nx.Graph, nx.DiGraph], nodes: tp.List, weight="weight"
) -> np.ndarray:
    weight_fn = _weight_function(graph, weight)
    weight_arr = np.empty((len(nodes), len(nodes)), dtype=np.uint16)
    weight_arr[:] = np.iinfo(np.uint16).max
    for edge in graph.edges:
        nd1_idx = nodes.index(edge[0])
        nd2_idx = nodes.index(edge[1])
        w = weight_fn(edge[0], edge[1], graph.edges[edge[0], edge[1]])
        weight_arr[nd1_idx, nd2_idx] = w
        weight_arr[nd2_idx, nd1_idx] = w
    return weight_arr


def _nxcolortonp(
    color: tp.Callable[[tp.Any], tp.FrozenSet[int]], nodes: tp.List
) -> np.ndarray:
    colors = np.unique(list(reduce(set.union, [set(color(node)) for node in nodes])))
    color_arr = np.zeros((len(nodes), colors.shape[0]), dtype=np.ubyte)
    for node in nodes:
        n_index = nodes.index(node)
        for node_color in color(node):
            color_arr[n_index, np.flatnonzero(colors == node_color).item()] = True
    return color_arr


def _nxgraphtomap(weight_arr: np.ndarray) -> tp.Mapping:
    cgraph = defaultdict(list)
    rows, cols = np.nonzero(weight_arr < np.iinfo(np.uint16).max)
    for i, j in zip(rows, cols):
        cgraph[int(i)].append(int(j))

    return colorful_algorithm.graphdict2map(dict(cgraph))


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
