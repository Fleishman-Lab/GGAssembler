from collections import defaultdict
from functools import partial
from heapq import heappop, heappush
from itertools import count, repeat

import networkx as nx
from networkx.algorithms.shortest_paths.weighted import _weight_function


def nested_defaultdict(default_factory, depth=1):
    result = partial(defaultdict, default_factory)
    for _ in repeat(None, depth - 1):
        result = partial(defaultdict, result)
    return result()


def all_shortest_paths(
    G, source, target, len_cutoff, weight=None, method="dijkstra", cutoff=None
):
    """Compute all shortest paths in the graph.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path.

    target : node
       Ending node for path.

    len_cutoff: int
        The maximum path length

    weight : None or string, optional (default = None)
       If None, every edge has weight/distance/cost 1.
       If a string, use this edge attribute as the edge weight.
       Any edge attribute not present defaults to 1.

    method : string, optional (default = 'dijkstra')
       The algorithm to use to compute the path lengths.
       Supported options: 'dijkstra'.
       Other inputs produce a ValueError.

    cutoff: numeric, optional (default = None)
        Maximum path weight

    Returns
    -------
    paths : generator of lists
        A generator of all paths between source and target.

    Raises
    ------
    ValueError
        If `method` is not among the supported options.

    NetworkXNoPath
        If `target` cannot be reached from `source`.


    Notes
    -----
    There may be many shortest paths between the source and target.

    """
    method = "unweighted" if weight is None else method
    if method == "dijkstra":
        pred, dist = dijkstra_predecessor_and_distance(
            G, source, weight=weight, cutoff=cutoff, len_cutoff=len_cutoff
        )
    else:
        raise ValueError("method not supported: {}".format(method))

    if target not in pred:
        raise nx.NetworkXNoPath(
            "Target {} cannot be reached" "from Source {}".format(target, source)
        )

    for path_len in pred[target].keys():
        stack = [[target, 0, path_len]]
        top = 0
        while top >= 0:
            node, i, l = stack[top]
            if node == source:
                yield (
                    [p for p, n, l in reversed(stack[: top + 1])],
                    dist[target][path_len],
                )
            if len(pred[node][l]) > i:
                top += 1
                if top == len(stack):
                    stack.append([pred[node][l][i], 0, l - 1])
                else:
                    stack[top] = [pred[node][l][i], 0, l - 1]
            else:
                stack[top - 1][1] += 1
                top -= 1


def dijkstra_predecessor_and_distance(
    G, source, len_cutoff=None, cutoff=None, weight="weight"
):
    """Compute weighted shortest path length and predecessors.

    Uses Dijkstra's Method to obtain the shortest weighted paths
    and return dictionaries of predecessors for each node and
    distance for each node from the `source`.

    Parameters
    ----------
    G : NetworkX graph

    source : node label
       Starting node for path

    len_cutoff: int
        The maximum path length

    cutoff : integer or float, optional
       Depth to stop the search. Only return paths with length <= cutoff.

    weight : string or function
       If this is a string, then edge weights will be accessed via the
       edge attribute with this key (that is, the weight of the edge
       joining `u` to `v` will be ``G.edges[u, v][weight]``). If no
       such edge attribute exists, the weight of the edge is assumed to
       be one.

       If this is a function, the weight of an edge is the value
       returned by the function. The function must accept exactly three
       positional arguments: the two endpoints of an edge and the
       dictionary of edge attributes for that edge. The function must
       return a number.

    Returns
    -------
    pred, distance : dictionaries
       Returns two dictionaries representing a list of predecessors
       of a node and the distance to each node.
       Warning: If target is specified, the dicts are incomplete as they
       only contain information for the nodes along a path to target.

    Raises
    ------
    NodeNotFound
        If `source` is not in `G`.

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    The list of predecessors contains more than one element only when
    there are more than one shortest paths to the key node.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.path_graph(5, create_using = nx.DiGraph())
    >>> pred, dist = nx.dijkstra_predecessor_and_distance(G, 0)
    >>> sorted(pred.items())
    [(0, []), (1, [0]), (2, [1]), (3, [2]), (4, [3])]
    >>> sorted(dist.items())
    [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]

    >>> pred, dist = nx.dijkstra_predecessor_and_distance(G, 0, 1)
    >>> sorted(pred.items())
    [(0, []), (1, [0])]
    >>> sorted(dist.items())
    [(0, 0), (1, 1)]
    """

    weight = _weight_function(G, weight)
    pred = nested_defaultdict(list, 2)  # dictionary of predecessors
    pred[source][1] = [source]
    return (
        pred,
        _dijkstra(G, source, weight, pred=pred, cutoff=cutoff, len_cutoff=len_cutoff),
    )


def _dijkstra(
    G, source, weight, len_cutoff, pred=None, paths=None, cutoff=None, target=None
):
    """Uses Dijkstra's algorithm to find shortest weighted paths from a
    single source.

    This is a convenience function for :func:`_dijkstra_multisource`
    with all the arguments the same, except the keyword argument
    `sources` set to ``[source]``.

    """
    return _dijkstra_multisource(
        G,
        [source],
        weight,
        pred=pred,
        paths=paths,
        cutoff=cutoff,
        target=target,
        len_cutoff=len_cutoff,
    )


def _dijkstra_multisource(
    G, sources, weight, len_cutoff, pred=None, paths=None, cutoff=None, target=None
):
    """Uses Dijkstra's algorithm to find shortest weighted paths

    Parameters
    ----------
    G : NetworkX graph

    sources : non-empty iterable of nodes
        Starting nodes for paths. If this is just an iterable containing
        a single node, then all paths computed by this function will
        start from that node. If there are two or more nodes in this
        iterable, the computed paths may begin from any one of the start
        nodes.

    weight: function
        Function with (u, v, data) input that returns that edges weight

    pred: dict of lists, optional(default=None)
        dict to store a list of predecessors keyed by that node
        If None, predecessors are not stored.

    paths: dict, optional (default=None)
        dict to store the path list from source to each node, keyed by node.
        If None, paths are not stored.

    target : node label, optional
        Ending node for path. Search is halted when target is found.

    cutoff : integer or float, optional
        Depth to stop the search. Only return paths with length <= cutoff.

    Returns
    -------
    distance : dictionary
        A mapping from node to shortest distance to that node from one
        of the source nodes.

    Raises
    ------
    NodeNotFound
        If any of `sources` is not in `G`.

    Notes
    -----
    The optional predecessor and path dictionaries can be accessed by
    the caller through the original pred and paths objects passed
    as arguments. No need to explicitly return pred or paths.

    """
    G_succ = G._succ if G.is_directed() else G._adj

    push = heappush
    pop = heappop
    dist = defaultdict(dict)  # dictionary of final distances
    seen = defaultdict(dict)
    # fringe is heapq with 4-tuples (distance, length, c, node)
    # use the count c to avoid comparing nodes (may not be able to)
    c = count()
    fringe = []
    for source in sources:
        if source not in G:
            raise nx.NodeNotFound("Source {} not in G".format(source))
        seen[source][0] = 0
        push(fringe, (0, 0, next(c), source))
    while fringe:
        (d, l, _, v) = pop(fringe)
        if v in dist and l in dist[v]:
            continue  # already searched this node.
        dist[v][l] = d
        if v == target:
            break
        for u, e in G_succ[v].items():
            cost = weight(v, u, e)
            if cost is None:
                continue
            vu_dist = dist[v][l] + cost
            if cutoff is not None:
                if vu_dist > cutoff:
                    continue
            vu_len = l + 1
            if len_cutoff is not None:
                if vu_len >= len_cutoff:
                    continue
            if u in dist and vu_len in dist[u]:
                if vu_dist < dist[u][vu_len]:
                    raise ValueError("Contradictory paths found:", "negative weights?")
            elif u not in seen or vu_len not in seen[u] or vu_dist < seen[u][vu_len]:
                seen[u][vu_len] = vu_dist
                push(fringe, (vu_dist, vu_len, next(c), u))
                if paths is not None:
                    paths[u][vu_len] = paths[v][l] + [u]
                if pred is not None:
                    pred[u][vu_len] = [v]
            elif vu_dist == seen[u][vu_len]:
                if pred is not None:
                    pred[u][vu_len].append(v)

    # The optional predecessor and path dictionaries can be accessed
    # by the caller via the pred and paths objects passed as arguments.
    return dist
