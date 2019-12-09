from itertools import permutations
from typing import List, Tuple, Dict, Callable, Union

import networkx as nx
import numpy as np

from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gate import Gate


SOURCE_NODE = "src"
SINK_NODE = "snk"


class GraphMaker:
    def __init__(self, gg_data):
        self.gg_data: GGData = gg_data

    def create_default_valid_edge_func(
        self,
        dna_var_poss: List[int],
        min_oligo_length: int,
        max_oligo_length: int,
        min_const_oligo_length: int,
        max_gate_crosstalk: int,
    ) -> Callable[[Gate, Gate], bool]:
        ggdata = self.gg_data
        dna_var_arr = np.array(dna_var_poss)

        def _is_valid_edge(nd1: Gate, nd2: Gate) -> bool:
            if nd1.idx + 3 >= nd2.idx:
                return False
            # segments with variable positions between them
            # are required to be at a certain length
            if np.any(np.logical_and(nd1.idx < dna_var_arr, dna_var_arr < nd2.idx)):
                if not (min_oligo_length < nd2.idx - nd1.idx < max_oligo_length):
                    return False
            else:
                if nd2.idx - nd1.idx < min_const_oligo_length:
                    return False
            if ggdata.gates_all_scores(nd1.bps, nd2.bps) > max_gate_crosstalk:
                return False
            return True

        return _is_valid_edge


def make_nodes(d_graph, dna: str, is_valid_node: Callable[[str, int], bool]) -> None:
    for ind in range(len(dna) - 3):
        gate_idxs = slice(ind, ind + 4)
        fcw = dna[gate_idxs]
        if is_valid_node(fcw, ind+1):
            d_graph.add_node(Gate(ind+1, fcw))


def create_default_valid_node_function(
    acceptable_fcws: List[str], var_dna_poss: List[int]
) -> Callable[[str, int], bool]:
    var_dna_arr = np.array(var_dna_poss)

    def _is_valid_node(fcw: str, ind: int) -> bool:
        if fcw not in acceptable_fcws:
            return False
        if np.any(np.logical_and(ind <= var_dna_arr, var_dna_arr <= ind + 3)):
            return False
        return True

    return _is_valid_node


def _add_source_sink(
    d_graph: nx.DiGraph, dna: str, var_poss: List[int]
) -> Tuple[Gate, Gate]:
    src = Gate(-1, SOURCE_NODE)
    snk = Gate(len(dna), SINK_NODE)
    d_graph.add_node(src)
    d_graph.add_node(snk)
    for nd1 in d_graph.nodes:
        if nd1.idx < var_poss[0]:
            d_graph.add_edge(src, nd1, weight=0)
        elif nd1.idx > var_poss[-1]:
            d_graph.add_edge(nd1, snk, weight=0)
    return src, snk


def create_default_weight_func(
    dna_pos_n_codons: Dict[int, List[str]]
) -> Callable[[Gate, Gate], int]:
    var_pos_arr = np.array(list(dna_pos_n_codons.keys()))

    def edge_weight(nd1: Gate, nd2: Gate) -> int:
        base_cost = 0
        if np.any(np.logical_and(nd1.idx < var_pos_arr, var_pos_arr < nd2.idx)):
            base_cost = 1
        return 1 + (nd2.idx - nd1.idx) * np.product(
            [
                len(codons) if nd1.idx < pos < nd2.idx else base_cost
                for pos, codons in dna_pos_n_codons.items()
            ]
        )

    return edge_weight


def make_edges(
    d_graph: nx.DiGraph,
    is_valid_edge: Callable[[Gate, Gate], bool],
    edge_weight: Callable[[Gate, Gate], Union[float, int]],
) -> None:
    for nd1, nd2 in permutations(d_graph.nodes, 2):
        if is_valid_edge(nd1, nd2):
            d_graph.add_edge(nd1, nd2, weight=edge_weight(nd1, nd2))


def build_custom_graph(
    gm: GraphMaker,
    dna: str,
    var_poss: List[int],
    is_valid_node: Callable[[str, int], bool],
    is_valid_edge: Callable[[Gate, Gate], bool],
    edge_weight: Callable[[Gate, Gate], Union[float, int]],
) -> Tuple[nx.DiGraph, Gate, Gate]:
    d_graph: nx.DiGraph = nx.DiGraph()
    make_nodes(d_graph, dna, is_valid_node)
    make_edges(d_graph, is_valid_edge, edge_weight)
    src, snk = _add_source_sink(d_graph, dna, var_poss)
    return d_graph, src, snk


def make_default_graph(
    gm: GraphMaker,
    dna: str,
    var_poss: List[int],
    dna_pos_n_codons: Dict[int, List[str]],
    min_oligo_length: int,
    max_oligo_length: int,
    min_const_oligo_length: int,
    gate_self_binding_min: int = 2000,
    gate_crosstalk_max: int = 1000,
) -> Tuple[nx.DiGraph, Gate, Gate]:
    acceptable_fcws = gm.gg_data.filter_self_binding_gates(gate_self_binding_min)
    is_valid_node = create_default_valid_node_function(acceptable_fcws, var_poss)
    is_valid_edge = gm.create_default_valid_edge_func(
        var_poss,
        min_oligo_length,
        max_oligo_length,
        min_const_oligo_length,
        gate_crosstalk_max,
    )
    edge_weight = create_default_weight_func(dna_pos_n_codons)
    return build_custom_graph(
        gm, dna, var_poss, is_valid_node, is_valid_edge, edge_weight
    )
