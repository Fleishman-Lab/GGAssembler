from itertools import permutations
from typing import List, Tuple

import networkx as nx

from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gate import Gate


SOURCE_NODE = "src"
SINK_NODE = "snk"

class GraphMaker:
    def __init__(self, gg_data):
        self.gg_data: GGData = gg_data

    def make_grpah(
        self,
        dna: str,
        const_poss: List[int],
        var_poss: List[int],
        min_oligo_length: int,
        max_oligo_length: int,
        min_const_oligo_length: int,
    ) -> nx.DiGraph:
        d_graph: nx.DiGraph = nx.DiGraph()
        self.make_nodes(d_graph, dna, const_poss)
        self.make_edges(
            d_graph,
            var_poss,
            min_oligo_length,
            max_oligo_length,
            min_const_oligo_length,
        )
        self.add_source_sink(
            d_graph,
            dna,
            var_poss,
            min_oligo_length,
            max_oligo_length,
            min_const_oligo_length,
        )
        return d_graph

    def make_nodes(self, d_graph, dna: str, const_dna_poss: List[int]) -> None:
        acceptable_fcws = self.gg_data.filter_self_binding_gates()
        for ind in range(len(dna) - 3):
            fcw = dna[ind : ind + 4]
            if _is_node(fcw, ind, acceptable_fcws, const_dna_poss):
                d_graph.add_node(Gate(ind, fcw))

    def make_edges(
        self,
        d_graph: nx.DiGraph,
        var_poss: List[int],
        min_oligo_length,
        max_oligo_length,
        min_const_oligo_length,
    ) -> None:
        for nd1, nd2 in permutations(d_graph.nodes, 2):
            if self._is_edge(
                nd1,
                nd2,
                var_poss,
                min_oligo_length,
                max_oligo_length,
                min_const_oligo_length,
            ):
                d_graph.add_edge(nd1, nd2, weight=nd2.index - nd1.index + 3)

    def add_source_sink(
        self,
        d_graph: nx.DiGraph,
        dna: str,
        var_poss: List[int],
        min_oligo_length,
        max_oligo_length,
        min_const_oligo_length,
    ) -> None:
        src = Gate(-1, SOURCE_NODE)
        snk = Gate(len(dna), SINK_NODE)
        d_graph.add_node(src)
        d_graph.add_node(snk)
        for nd1 in d_graph.nodes:
            # if theres a variable position before nd1, conncet it so source only if it is short enough
            # otherwise connect it anyway (constant segment)
            if any([-1 < p < nd1.index for p in var_poss]):
                if min_oligo_length < nd1.index < max_oligo_length:
                    d_graph.add_edge(src, nd1)
            else:
                d_graph.add_edge(src, nd1)
            # if theres a variable position between nd1 and sink, conncet it so sink only if it is short enough
            # otherwise connect it anyway (constant segment)
            if any([nd1.index < p for p in var_poss]):
                if min_oligo_length < len(dna) - nd1.index < max_oligo_length:
                    d_graph.add_edge(nd1, snk)
            else:
                d_graph.add_edge(nd1, snk)

    def _is_edge(
        self,
        nd1: Gate,
        nd2: Gate,
        dna_var_poss: List[int],
        min_oligo_length,
        max_oligo_length,
        min_const_oligo_length,
    ) -> bool:
        if nd1.index + 3 >= nd2.index:
            return False

        # segments with variable positions between them are required to be at a certain length
        if any([nd1.index < p < nd2.index for p in dna_var_poss]):
            if not (min_oligo_length < nd2.index - nd1.index + 3 < max_oligo_length):
                return False
        else:
            if nd2.index - nd1.index + 3 < min_const_oligo_length:
                return False
        if self.gg_data.gates_all_scores(nd1.bps, nd2.bps) > 1000:
            return False
        return True


def _is_node(
    fcw: str, ind: int, acceptable_fcws: List[str], const_dna_poss: List[int]
) -> bool:
    if fcw not in acceptable_fcws:
        return False
    if any([p not in const_dna_poss for p in range(ind, ind + 4)]):
        return False
    return True
