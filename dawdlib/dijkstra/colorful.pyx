# distutils: language = c++
# distutils: extra_compile_args=-Ofast -std=c++11
# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=False, overflowcheck.fold=False

# '''
# # cython: infer_types=True
# '''
# '''
# # cython: binding=True
# # cython: profile=True
# # cython: linetrace=True
# # distutils: define_macros=CYTHON_TRACE_NOGIL=1
# '''
import numpy as np
import typing as tp
from itertools import chain
import networkx as nx

# from dawdlib.golden_gate.gate_data import GGData
# from dawdlib.dijkstra.colorful_utils import color_gates, nxtonumpy, return_shortest_path, yield_shortest_paths

from libcpp.map cimport map
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t
from libcpp.limits cimport numeric_limits
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free, rand
from libc.string cimport memset

ctypedef fused umy_type:
    uint8_t
    uint16_t
    uint32_t
    uint64_t

ctypedef fused my_type:
    int8_t
    int16_t
    int32_t
    int64_t

cdef extern from "cpp_priority_queue.hpp" nogil:
    cdef cppclass greater_priority_queue[T]:
        greater_priority_queue() except +
        greater_priority_queue(greater_priority_queue&) except +
        bint empty()
        void pop()
        void push(T&)
        size_t size()
        T& top()
        # C++11 methods
        void swap(greater_priority_queue&)


# @cython.embedsignature(True)
# def all_shortest_paths(ggdata, G, source, target, no_colors=0, len_cutoff=0, weight=None, method="dijkstra", retries = 1, res_graph: nx.Graph = None):
#     """Compute all shortest paths in the graph.
#
#     Parameters
#     ----------
#     ggdata: GGData
#        GGdata for the reaction.
#
#     G : NetworkX graph
#
#     source : node
#        Starting node for path.
#
#     target : node
#        Ending node for path.
#
#     no_colors: int
#         The number of colors to use in the random recoloring.
#         The currently supported maximum is set to 18.
#
#     len_cutoff: int, optional (default = 0)
#         The maximum path length (numbe of gates in the solution)
#
#     weight : None or string, optional (default = None)
#        If None, every edge has weight/distance/cost 1.
#        If a string, use this edge attribute as the edge weight.
#        Any edge attribute not present defaults to 1.
#
#     method : string, optional (default = 'dijkstra')
#        The algorithm to use to compute the path lengths.
#        Supported options: 'dijkstra'.
#        Other inputs produce a ValueError.
#
#     retries: int, optional (default = 1)
#         The number of random coloring retries to check before raising NetworkXNoPath exception
#
#     res_graph: nx.Graph, optional (default = None)
#         Allows the user to provide a custom restriction graph used to color nodes
#
#     Returns
#     -------
#     paths : generator of lists
#         A generator of all paths between source and target.
#
#     Raises
#     ------
#     ValueError
#         If `method` is not among the supported options.
#         If both `no_colors` and `len_cutoff` are not provided.
#
#     NetworkXNoPath
#         If `target` cannot be reached from `source`.
#
#
#     Notes
#     -----
#     There may be many shortest paths between the source and target.
#
#     """
#     nodes, cgraph, weight_arr, sources, trgt, color_arr, new_color, no_colors, pred, s_pred, dist, seen = _prep_data(ggdata, G, source, target, no_colors=no_colors, len_cutoff=len_cutoff, weight=weight, method=method, retries = retries, res_graph = res_graph)
#     try:
#         return all_colorful_shortest_paths(
#             nodes, cgraph, weight_arr, sources, trgt, color_arr, new_color, no_colors, pred, s_pred, dist, seen, limit = len_cutoff
#         )
#     except nx.NetworkXNoPath:
#         raise nx.NetworkXNoPath(
#             "Target {} cannot be reached" "from Source {}".format(target, source)
#         )
#
#
# @cython.embedsignature(True)
# def shortest_path(ggdata, G, source, target, no_colors=0, len_cutoff=0, weight=None, method="dijkstra", retries = 1, res_graph: nx.Graph = None):
#     """Compute all shortest paths in the graph.
#
#     Parameters
#     ----------
#     ggdata: GGData
#        GGdata for the reaction.
#
#     G : NetworkX graph
#
#     source : node
#        Starting node for path.
#
#     target : node
#        Ending node for path.
#
#     no_colors: int
#         The number of colors to use in the random recoloring.
#         The currently supported maximum is set to 18.
#
#     len_cutoff: int, optional (default = 0)
#         The maximum path length
#
#     weight : None or string, optional (default = None)
#        If None, every edge has weight/distance/cost 1.
#        If a string, use this edge attribute as the edge weight.
#        Any edge attribute not present defaults to 1.
#
#     method : string, optional (default = 'dijkstra')
#        The algorithm to use to compute the path lengths.
#        Supported options: 'dijkstra'.
#        Other inputs produce a ValueError.
#
#     retries: int, optional (default = 1)
#         The number of random coloring retries to check before raising NetworkXNoPath exception
#
#     res_graph: nx.Graph, optional (default = None)
#         Allows the user to provide a custom restriction graph used to color nodes
#
#     Returns
#     -------
#     paths : generator of lists
#         A generator of all paths between source and target.
#
#     Raises
#     ------
#     ValueError
#         If `method` is not among the supported options.
#         If both `no_colors` and `len_cutoff` are not provided.
#
#     NetworkXNoPath
#         If `target` cannot be reached from `source`.
#
#
#     Notes
#     -----
#     There may be many shortest paths between the source and target.
#
#     """
#     nodes, cgraph, weight_arr, sources, trgt, color_arr, new_color, no_colors, pred, s_pred, dist, seen = _prep_data(ggdata, G, source, target, no_colors=no_colors, len_cutoff=len_cutoff, weight=weight, method=method, retries = retries, res_graph = res_graph)
#     try:
#         return colorful_shortest_path(
#             nodes, cgraph, weight_arr, sources, trgt, color_arr, new_color, no_colors, pred, s_pred, dist, seen, limit = len_cutoff
#         )
#     except nx.NetworkXNoPath:
#         raise nx.NetworkXNoPath(
#             "Target {} cannot be reached" "from Source {}".format(target, source)
#         )
#
#
# def _prep_data(ggdata, G, source, target, no_colors=0, len_cutoff=0, weight=None, method="dijkstra", retries = 1, res_graph: nx.Graph = None):
#     if method != 'dijkstra':
#         raise ValueError("method not supported: {}".format(method))
#     if not no_colors and not len_cutoff:
#         raise ValueError("Either no_colors or len_cutoff are required.")
#
#     gate_colors, num_of_colors = color_gates(G, ggdata, r_graph=res_graph)
#     if not no_colors:
#         no_colors = sum(num_of_colors[:len_cutoff])
#     if no_colors > 18:
#         print(f"Warning! number of colors required {no_colors} is larger than recommended.")
#     return nxtonumpy(G, [source], target, gate_colors.get, no_colors)
#
#
# @cython.embedsignature(True)
# def all_colorful_shortest_paths(
#     nodes: tp.List[int, tp.Any],
#     G: tp.Dict[int, tp.FrozenSet[int]],
#     weight: np.ndarray,
#     sources: np.ndarray,
#     target: int,
#     colors: np.ndarray,
#     color_map: np.ndarray,
#     no_colors: int,
#     pred: np.ndarray,
#     s_pred: np.ndarray,
#     dist: np.ndarray,
#     seen: np.ndarray,
#     limit: int = 0,
#     repetitions: int = 1,
# ) -> tp.Generator[tp.Any]:
#     """
#     :param list nodes: The list of nodes, used to identify each node with an index
#     :param dict G: A mapping representation of the golden gate graph
#     :param numpy.ndarray weight: A numpy array representations of graph edge weights
#     :param numpy.ndarray sources: A numpy array holding indices of sources
#     :param int target: The path target
#     :param numpy.ndarray colors: A numpy array mapping between a node and it's colors
#     :param numpy.ndarray color_map: A numpy array holding the random recoloring
#     :param int no_colors: How many colors were used to created the dist and seen arrays
#     :param numpy.ndarray pred: A numpy array pointing to the predecessor of each node in the paths found
#     :param numpy.ndarray s_pred: A numpy array pointing to a single predecessor of each node
#     :param numpy.ndarray dist: A numpy array holding the minimum distance from source to node using specific colors
#     :param numpy.ndarray seen: A numpy array holding the minimum distance from source to node using specific colors
#     :param int limit: A limit to the number of gates found in the shortest path
#     :param int repetitions: How many times to try and find a shortest path using a random recoloring of no_colors
#     :return: A generator yielding shortest colorful paths
#     :rtype: generator
#     """
#     found, res = find_shortest_paths(G, weight, sources, colors, color_map, no_colors, pred, s_pred, dist, seen, limit, repetitions)
#     if found:
#         return chain.from_iterable((yield_shortest_paths(nodes, src, target, color_map, pred) for src in sources))
#     else:
#         raise nx.NetworkXNoPath()
#
#
# @cython.embedsignature(True)
# def colorful_shortest_path(
#     nodes: tp.List[int, tp.Any],
#     G: tp.Dict[int, tp.FrozenSet[int]],
#     weight: np.ndarray,
#     sources: np.ndarray,
#     target: int,
#     colors: np.ndarray,
#     color_map: np.ndarray,
#     no_colors: int,
#     pred: np.ndarray,
#     s_pred: np.ndarray,
#     dist: np.ndarray,
#     seen: np.ndarray,
#     limit: int = 0,
#     repetitions: int = 1,
# ) -> tp.List[tp.Any]:
#     """
#     :param list nodes: The list of nodes, used to identify each node with an index
#     :param dict G: A mapping representation of the golden gate graph
#     :param numpy.ndarray weight: A numpy array representations of graph edge weights
#     :param numpy.ndarray sources: A numpy array holding indices of sources
#     :param int target: The path target
#     :param numpy.ndarray colors: A numpy array mapping between a node and it's colors
#     :param numpy.ndarray color_map: A numpy array holding the random recoloring
#     :param int no_colors: How many colors were used to created the dist and seen arrays
#     :param numpy.ndarray pred: A numpy array pointing to the predecessor of each node in the paths found
#     :param numpy.ndarray s_pred: A numpy array pointing to a single predecessor of each node
#     :param numpy.ndarray dist: A numpy array holding the minimum distance from source to node using specific colors
#     :param numpy.ndarray seen: A numpy array holding the minimum distance from source to node using specific colors
#     :param int limit: A limit to the number of gates found in the shortest path
#     :param int repetitions: How many times to try and find a shortest path using a random recoloring of no_colors
#     :return: A generator yielding shortest colorful paths
#     :rtype: generator
#     """
#     found, res = find_shortest_paths(G, weight, sources, colors, color_map, no_colors, pred, s_pred, dist, seen, limit, repetitions)
#     if found:
#         for src in sources:
#             yield return_shortest_path(nodes, src, target, s_pred, seen, color_map, limit)
#     else:
#         raise nx.NetworkXNoPath()


def find_shortest_paths(
    map[int, vector[int]] G,
    uint16_t[:, ::1] weight,
    int[::1] sources,
    uint8_t[:, ::1] colors,
    umy_type[::1] color_map,
    int no_colors,
    uint8_t[:, ::1] pred,
    umy_type[:, ::1] s_pred,
    umy_type[:, ::1] dist,
    umy_type[:, ::1] seen,
    int limit = 0,
    int repetitions = 1,
):
    cdef:
        int j = 0
        int i = 0
        int res = 0

    for i in range(repetitions):
        pred[:] = 0
        s_pred[:] = numeric_limits[umy_type].max()
        res = random_colorful_shortest_paths(G, weight, sources, colors, color_map, no_colors, pred, s_pred, dist, seen, limit)
        for j in range(pred.shape[1]):
            if pred[pred.shape[0]-1, j] > 0:
                return True, res
    return False, res



# def yield_colorful_shortest_paths(src, trgt, s_pred, seen, int limit=0):
#     seen_count = np.sum(seen < np.iinfo(seen.dtype).max, axis=1)
#     seen_argsort = np.argsort(seen, axis=1)
#
#     stack = [[trgt, 0, 0]]
#     top = 0
#     visited = defaultdict(list)
#
#     while top >= 0:
#         v, i, p_len = stack[top]
#         if v == src:
#             yield [p for p, _, _ in reversed(stack[: top + 1])]
#         if seen_count[v] > i and not 0 < limit < p_len:
#             # pred is:
#             u = s_pred[v, seen_argsort[v, i]]
#             if u in visited[(v, p_len)]:
#                 stack[top][1] += 1
#                 continue
#             visited[(v, p_len)].append(u)
#             top += 1
#             u_tup = [u, 0, p_len + 1]
#             if top == len(stack):
#                 stack.append(u_tup)
#             else:
#                 stack[top] = u_tup
#         else:
#             stack[top - 1][1] += 1
#             top -= 1


# def yield_colorful_shortest_paths(nodes, src, trgt, node_colors, pred, limit=0):
#     pred_dict = defaultdict(list)
#     for key, val in np.argwhere(pred).tolist():
#         pred_dict[key].append(val)
#     stack = [[trgt, 0, node_colors[trgt]]]
#     top = 0
#     while top >= 0:
#         node, i, color, p_len = stack[top]
#         if node == src:
#             yield [nodes[p] for p, n, c, l in reversed(stack[: top + 1])]
#         if len(pred_dict[node]) > i:
#             if not (color & node_colors[pred_dict[node][i]]) and not 0 < limit < p_len:
#                 top += 1
#                 if top == len(stack):
#                     stack.append(
#                         [
#                             pred_dict[node][i],
#                             0,
#                             color | node_colors[pred_dict[node][i]],
#                             p_len + 1,
#                         ]
#                     )
#                 else:
#                     stack[top] = [
#                         pred_dict[node][i],
#                         0,
#                         color | node_colors[pred_dict[node][i]],
#                         p_len + 1,
#                     ]
#             else:
#                 stack[top][1] += 1
#         else:
#             stack[top - 1][1] += 1
#             top -= 1


cdef int random_colorful_shortest_paths(
    map[int, vector[int]] G,
    uint16_t[:, ::1] weight,
    int[::1] sources,
    uint8_t[:, ::1] colors,
    umy_type[::1] color_map,
    int no_colors,
    uint8_t[:, ::1] pred,
    umy_type[:, ::1] s_pred,
    umy_type[:, ::1] dist,
    umy_type[:, ::1] seen,
    int limit = 0,
):
    cdef:
        int res = 0
    random_recolor(colors, color_map, no_colors)

    # Check if we found paths and return them as iterator
    res = predecessor_and_distance(G, weight, sources, color_map, pred, s_pred, dist, seen, limit)
    return res


cdef void random_recolor(uint8_t[:, ::1] colors, umy_type[::1] color_map, int32_t no_colors):
    cdef:
        uint8_t *random_coloring = <uint8_t *> malloc(colors.shape[0]*no_colors*sizeof(uint8_t))

    if not random_coloring:
        raise MemoryError()

    try:
        memset(random_coloring, 0, colors.shape[0]*no_colors*sizeof(uint8_t))
        recolor_c(colors, random_coloring, no_colors)
        boolmat2int_c(random_coloring, color_map,  no_colors)
    finally:
        free(random_coloring)

# @cython.embedsignature(True)
# def recolor(color_arr, no_colors):
#     cdef Py_ssize_t cshape = color_arr.shape[1]
#     recolor_arr = np.zeros((cshape, no_colors), np.ubyte)
#     rand_colors = np.random.randint(0, no_colors, color_arr.shape[1])
#     recolor_arr[np.arange(cshape), rand_colors] = True
#     return np.matmul(color_arr, recolor_arr)


@cython.cdivision(True)
cdef void recolor_c(uint8_t[:, ::1] color_arr, uint8_t  *new_color_arr, int32_t no_colors):
    cdef:
        int current_colors = color_arr.shape[1]
        int no_nodes = color_arr.shape[0]
        int i, j, color
        uint8_t *recolor_arr = <uint8_t *> malloc(current_colors*sizeof(uint8_t))

    if not recolor_arr:
        raise MemoryError()

    try:
        memset(recolor_arr, 0, current_colors*sizeof(uint8_t))
        for i in range(current_colors):
            recolor_arr[i] = rand() % no_colors
        for i in range(no_nodes):
            for j in range(current_colors):
                if color_arr[i , j] == 1:
                    new_color_arr[i * no_colors + recolor_arr[j]] = 1
    finally:
        free(recolor_arr)

# @cython.embedsignature(True)
# cdef umy_type bool2int(uint8_t[:] arr, umy_type* res_ref) nogil:
#     cdef Py_ssize_t x_max = arr.shape[0]
#     cdef Py_ssize_t x
#     cdef umy_type res = deref(res_ref)
#     res = 0
#     for x in range(x_max):
#         res = (res << 1) ^ arr[x]
#     return res


# @cython.embedsignature(True)
# def boolmat2int(uint8_t[:, :] arr, umy_type[::1] result):
#     cdef:
#         Py_ssize_t x_max = arr.shape[0]
#         Py_ssize_t x
#
#     for x in range(x_max):
#         result[x] = bool2int(arr[x], &result[x])
#     # return result


cdef umy_type bool2int_c(uint8_t *arr, umy_type res, int32_t size) nogil:
    # cdef Py_ssize_t x_max = arr.shape[0]
    cdef int x
    res = 0
    for x in range(size):
        res = (res << 1) ^ arr[x]
    return res


cdef void boolmat2int_c(uint8_t *arr, umy_type[::1] result, int32_t no_colors) nogil:
    cdef:
        Py_ssize_t x_max = result.shape[0]
        Py_ssize_t x

    result[:] = 0
    for x in range(x_max):
        result[x] = bool2int_c(&arr[x*no_colors], result[x], no_colors)


@cython.embedsignature(True)
def graphdict2map(dict graph):
    cdef int v, u
    cdef map[int, vector[int]] *cgraph_ref = new map[int, vector[int]]()
    cdef map[int, vector[int]] cgraph = deref(cgraph_ref)
    cdef vector[int] *vec_ref
    cdef vector[int] vec
    cdef Py_ssize_t vec_size
    for v, val in graph.items():
        vec_size = len(val)
        vec_ref = new vector[int](vec_size)
        vec = deref(vec_ref)
        for u in range(vec_size):
            vec[u] = val[u]
        cgraph[v] = vec

    return cgraph


# @cython.embedsignature(True)
# cdef void pop_dist_seen(map[int, map[int, umy_type]] *d_ref,  map[int, map[int, umy_type]] *s_ref, int no_nodes):
#     cdef:
#         map[int, map[int, umy_type]] d = deref(d_ref)
#         map[int, map[int, umy_type]] s = deref(s_ref)
#         int i
#     for i in range(no_nodes):
#         d[i] = map[int, umy_type]()
#         s[i] = map[int, umy_type]()


# @cython.embedsignature(True)
# def colorsdict2map(dict color_dict):
#     cdef int v, val
#     cdef map[int, int] *color_ref = new map[int, int]()
#     cdef map[int, int] color = deref(color_ref)
#     cdef Py_ssize_t vec_size
#     for v, val in color_dict.items():
#         color[v] = val
#     return color


cdef vector[umy_type] vectorize(umy_type d, umy_type p_len, umy_type v_col, umy_type v):
    cdef:
        vector[umy_type] vec = vector[umy_type](4)

    vec[3] = v
    vec[2] = v_col
    vec[1] = p_len
    vec[0] = d
    return vec


# @cython.embedsignature(True)
cdef int _multisource(
    map[int, vector[int]] G,
    uint16_t[:, ::1] weight,
    int[::1] sources,
    umy_type[::1] color_map,
    uint8_t[:, ::1] pred,
    umy_type[:, ::1] s_pred,
    umy_type[:, ::1] dist,
    umy_type[:, ::1] seen,
    int limit=0
):
    max_sources = sources.shape[0]
    # fringe is heapq with 4-tuples (distance, length, color, node)
    cdef:
        umy_type vu_dist, d, source_color, source, v, v_col, vu_col, u_col, i, u, p_len
        int counter = 0
        vector[umy_type] *vec_ref
        vector[umy_type] vec
        greater_priority_queue[vector[umy_type]] fringe
    for i in range(max_sources):
        source = sources[i]
        source_color = color_map[source]
        seen[source, source_color] = 0
        vec = vectorize(0, 0, source_color, source)
        fringe.push(vec)
    while fringe.size() != 0:
        counter += 1
        vec = fringe.top()
        d = vec[0]
        p_len = vec[1]
        v_col = vec[2]
        v = vec[3]
        fringe.pop()
        if 0 < limit < p_len:
            continue
        if numeric_limits[umy_type].min() < dist[v, v_col]:
            continue
        dist[v, v_col] = d
        for u in G[v]:
            u_col = color_map[u]
            if v_col & u_col:
                continue
            vu_col = v_col | u_col
            vu_dist = d + weight[v, u]
            if vu_dist < dist[u, vu_col]:
                return -counter
            if vu_dist < seen[u, vu_col]:
                seen[u, vu_col] = vu_dist
                vec = vectorize(vu_dist, p_len + 1, vu_col, u)
                fringe.push(vec)
                pred[u, v] = 1
                s_pred[u, vu_col] = v
            elif vu_dist == seen[u, vu_col]:
                pred[u, v] = 1
                s_pred[u, vu_col] = v
    return counter


# @cython.embedsignature(True)
cdef int predecessor_and_distance(
    map[int, vector[int]] G,
    uint16_t[:, ::1] weight,
    int[::1] sources,
    umy_type[::1] color,
    uint8_t[:, ::1] pred,
    umy_type[:, ::1] s_pred,
    umy_type[:, ::1] dist,
    umy_type[:, ::1] seen,
    int limit=0
):
    cdef:
        int res

    dist[:] = numeric_limits[umy_type].min()
    seen[:] = numeric_limits[umy_type].max()
    res  = _multisource(
            G, weight, sources, color, pred, s_pred, dist, seen, limit
    )
    return res


# TODO: Extra code to add support for a really large number of colors ( > 25)
#  using mappings instead of numpy arrays.This comes at a high cost of computation time
#
#
# cdef int d_random_colorful_shortest_paths(
#     map[int, vector[int]] G,
#     my_type[:, ::1] weight,
#     int[::1] sources,
#     uint8_t[:, ::1] colors,
#     my_type[::1] color_map,
#     int no_colors,
#     uint8_t[:, ::1] pred,
#     int limit = 0,
# ):
#     cdef:
#         map[int, map[int, my_type]] *d_ref = new map[int, map[int, my_type]]()
#         map[int, map[int, my_type]] *s_ref = new map[int, map[int, my_type]]()
#         int res = 0
#
#     random_recolor(colors, color_map, no_colors)
#
#     # Check if we found paths and return them as iterator
#     res = d_predecessor_and_distance(G, weight, sources, color_map, pred, d_ref, s_ref, limit)
#     return res
#
#

# @cython.embedsignature(True)
# cdef int _d_multisource(
#     map[int, vector[int]] G,
#     my_type[:, ::1] weight,
#     int[::1] sources,
#     my_type[::1] color_map,
#     uint8_t[:, ::1] pred,
#     map[int, map[int, my_type]] *d_ref,
#     map[int, map[int, my_type]] *s_ref,
#     int limit=0
# ):
#     max_sources = sources.shape[0]
#     # fringe is heapq with 3-tuples (distance, color, node)
#     cdef:
#         my_type d, vu_dist, source, v, v_col, vu_col, u_col, source_color, i, u, p_len
#         int counter = 0
#         vector[my_type] *vec_ref
#         vector[my_type] vec
#         priority_queue[vector[my_type]] fringe
#         map[int, map[int, my_type]] dist = deref(d_ref)
#         map[int, map[int, my_type]] seen = deref(s_ref)
#     for i in range(max_sources):
#         source = sources[i]
#         source_color = color_map[source]
#         seen[source][source_color] = 0
#         vec = vectorize(0, 0, source_color, source)
#         fringe.push(vec)
#     while fringe.size() != 0:
#         counter += 1
#         vec = fringe.top()
#         d = -vec[0]
#         p_len = vec[1]
#         v_col = vec[2]
#         v = vec[3]
#         fringe.pop()
#         if limit > 0 and p_len > limit:
#             continue
#         if 0 < dist[v][v_col]:
#             continue
#         dist[v][v_col] = d
#         for u in G[v]:
#             u_col = color_map[u]
#             if v_col & u_col:
#                 continue
#             vu_col = v_col | u_col
#             vu_dist = d + weight[v, u]
#             if vu_dist < dist[u][vu_col]:
#                 return -counter
#             if 0 == seen[u][vu_col] or vu_dist < seen[u][vu_col]:
#                 seen[u][vu_col] = vu_dist
#                 vec = vectorize(-vu_dist, p_len + 1, vu_col, u)
#                 fringe.push(vec)
#                 pred[u, v] = 1
#             elif vu_dist == seen[u][vu_col]:
#                 pred[u, v] = 1
#     # The optional predecessor and path dictionaries can be accessed
#     # by the caller via the pred and paths objects passed as arguments.
#     return counter

# @cython.embedsignature(True)
# def d_predecessor_and_distance(
#     map[int, vector[int]] G,
#     my_type[:, ::1] weight,
#     int[::1] sources,
#     my_type[::1] color,
#     uint8_t[:, ::1] pred,
#     map[int, map[int, my_type]] *d_ref,
#     map[int, map[int, my_type]] *s_ref,
#     int limit=0
# ):
#     cdef:
#         int res
#
#     pop_dist_seen(d_ref, s_ref, weight.shape[0])
#     res  = _d_multisource(
#             G, weight, sources, color, pred, d_ref, s_ref, limit
#     )
#     return res
