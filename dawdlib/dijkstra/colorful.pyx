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

