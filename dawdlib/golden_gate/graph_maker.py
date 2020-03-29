from itertools import combinations
from typing import Callable, Dict, List, Tuple, Union

import networkx as nx
import numpy as np

from dawdlib.golden_gate.constants import (
    CONST_COST,
    SINK_NODE,
    SOURCE_NODE,
    VAR_ADD_COST,
)
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.utils import Requirements, syn_muts

SOURCE = Gate(None, SOURCE_NODE)
TARGET = Gate(None, SINK_NODE)


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
            if nd1.overlap(nd2):
                return False
            # segments with variable positions between them
            # are required to be at a certain length
            if np.any(np.logical_and(nd1.idx < dna_var_arr, dna_var_arr < nd2.idx)):
                # 4 was added to account for nd2 length (4) plus the fact we need to include the entire gate
                if not (min_oligo_length < nd2 - nd1 < max_oligo_length):
                    return False
            else:
                if nd2 - nd1 < min_const_oligo_length:
                    return False
            if ggdata.gates_all_scores(nd1.bps, nd2.bps) > max_gate_crosstalk:
                return False
            return True

        return _is_valid_edge


def make_nodes(
    d_graph, dna: str, is_valid_node: Callable[[Gate], bool], add_syn: bool = False
) -> None:
    if add_syn and len(dna) % 3 == 0:
        add_syn = True
    possible_nodes: List[Gate] = []
    for ind in range(len(dna) - 3):
        gate_idxs = slice(ind, ind + 4)
        fcw = dna[gate_idxs]
        possible_nodes.append(Gate(ind, fcw, False))
    if add_syn:
        possible_nodes += syn_muts(dna)
    d_graph.add_nodes_from(filter(is_valid_node, possible_nodes))


def create_default_valid_node_function(
    acceptable_fcws: List[str], var_dna_poss: List[int]
) -> Callable[[Gate], bool]:
    var_dna_arr = np.array(var_dna_poss)

    def _is_valid_node(gate: Gate) -> bool:

        if gate.bps not in acceptable_fcws:
            return False
        if np.any(np.logical_and(gate.idx <= var_dna_arr, var_dna_arr <= gate.idx + 3)):
            return False
        return True

    return _is_valid_node


def _add_source_sink(
    d_graph: nx.DiGraph, dna: str, var_poss: List[int]
) -> Tuple[Gate, Gate]:
    edges = []
    for nd1 in d_graph.nodes:
        if nd1.idx < var_poss[0]:
            edges.append((SOURCE, nd1, 0))
        elif nd1.idx > var_poss[-1]:
            edges.append((nd1, TARGET, 0))
    d_graph.add_node(SOURCE)
    d_graph.add_node(TARGET)
    d_graph.add_weighted_edges_from(edges)
    return SOURCE, TARGET


def create_default_weight_func(
    dna_pos_n_codons: Dict[int, List[str]],
    oligo_addition: int = VAR_ADD_COST,
    const_cost: int = CONST_COST,
) -> Callable[[Gate, Gate], int]:
    var_pos_arr = np.array(list(dna_pos_n_codons.keys()))

    def edge_weight(nd1: Gate, nd2: Gate) -> int:
        if np.any(np.logical_and(nd1.idx < var_pos_arr, var_pos_arr < nd2.idx)):
            return (nd2 - nd1 + oligo_addition) * np.product(
                [
                    len(codons) if nd1.idx < pos < nd2.idx else 1
                    for pos, codons in dna_pos_n_codons.items()
                ]
            )
        return const_cost

    return edge_weight


def make_edges(
    d_graph: nx.DiGraph,
    is_valid_edge: Callable[[Gate, Gate], bool],
    edge_weight: Callable[[Gate, Gate], Union[float, int]],
) -> None:
    d_graph.add_weighted_edges_from(
        map(
            lambda x: x + (edge_weight(*x),),
            filter(
                lambda x: is_valid_edge(*x),
                map(
                    lambda x: x if x[0].idx <= x[1].idx else (x[1], x[0]),
                    combinations(d_graph.nodes, 2),
                ),
            ),
        )
    )


def build_custom_graph(
    dna: str,
    var_poss: List[int],
    is_valid_node: Callable[[Gate], bool],
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
    reqs: Requirements,
) -> Tuple[nx.DiGraph, Gate, Gate]:
    acceptable_fcws = gm.gg_data.filter_self_binding_gates(reqs.gate_self_binding_min)
    is_valid_node = create_default_valid_node_function(acceptable_fcws, var_poss)
    is_valid_edge = gm.create_default_valid_edge_func(
        var_poss,
        reqs.min_oligo_length,
        reqs.max_oligo_length,
        reqs.min_const_oligo_length,
        reqs.gate_crosstalk_max,
    )
    edge_weight = create_default_weight_func(
        dna_pos_n_codons, reqs.oligo_addition, reqs.const_cost
    )
    return build_custom_graph(dna, var_poss, is_valid_node, is_valid_edge, edge_weight)
