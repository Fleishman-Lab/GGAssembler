from itertools import permutations
from typing import List, Tuple

import networkx as nx

from dawdlib.gate_data.gate_data import GGData
from dawdlib.gate_data.gate import Gate


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
        MIN_OLIGO_LENGTH,
        MAX_OLIGO_LENGTH,
        MIN_CONST_OLIGO_LENGTH,
    ) -> nx.DiGraph:
        d_graph: nx.DiGraph = nx.DiGraph()
        self.make_nodes(d_graph, dna, const_poss)
        self.make_edges(
            d_graph,
            var_poss,
            MIN_OLIGO_LENGTH,
            MAX_OLIGO_LENGTH,
            MIN_CONST_OLIGO_LENGTH,
        )
        self.add_source_sink(
            d_graph,
            dna,
            var_poss,
            MIN_OLIGO_LENGTH,
            MAX_OLIGO_LENGTH,
            MIN_CONST_OLIGO_LENGTH,
        )
        return d_graph

    def make_nodes(self, d_graph, dna: str, const_dna_poss: List[int]) -> None:
        acceptable_fcws = self.gg_data.get_all_self_binding_gates()
        for ind in range(len(dna) - 3):
            fcw = dna[ind : ind + 4]
            if _is_node(fcw, ind, acceptable_fcws, const_dna_poss):
                d_graph.add_node(Gate(ind, fcw))

    def make_edges(
        self,
        d_graph: nx.DiGraph,
        var_poss: List[int],
        MIN_OLIGO_LENGTH,
        MAX_OLIGO_LENGTH,
        MIN_CONST_OLIGO_LENGTH,
    ) -> None:
        for nd1, nd2 in permutations(d_graph.nodes, 2):
            if self._is_edge(
                nd1,
                nd2,
                var_poss,
                MIN_OLIGO_LENGTH,
                MAX_OLIGO_LENGTH,
                MIN_CONST_OLIGO_LENGTH,
            ):
                d_graph.add_edge(nd1, nd2)

    def add_source_sink(
        self,
        d_graph: nx.DiGraph,
        dna: str,
        var_poss: List[int],
        MIN_OLIGO_LENGTH,
        MAX_OLIGO_LENGTH,
        MIN_CONST_OLIGO_LENGTH,
    ) -> None:
        src = Gate(-1, SOURCE_NODE)
        snk = Gate(len(dna), SINK_NODE)
        d_graph.add_node(src)
        d_graph.add_node(snk)
        for nd1 in d_graph.nodes:
            # if theres a variable position before nd1, conncet it so source only if it is short enough
            # otherwise connect it anyway (constant segment)
            if any([-1 < p < nd1.index for p in var_poss]):
                if MIN_OLIGO_LENGTH < nd1.index < MAX_OLIGO_LENGTH:
                    d_graph.add_edge(src, nd1)
            else:
                d_graph.add_edge(src, nd1)
            # if theres a variable position between nd1 and sink, conncet it so sink only if it is short enough
            # otherwise connect it anyway (constant segment)
            if any([nd1.index < p for p in var_poss]):
                if MIN_OLIGO_LENGTH < len(dna) - nd1.index < MAX_OLIGO_LENGTH:
                    d_graph.add_edge(nd1, snk)
            else:
                d_graph.add_edge(nd1, snk)

    def _is_edge(
        self,
        nd1: Gate,
        nd2: Gate,
        dna_var_poss: List[int],
        MIN_OLIGO_LENGTH,
        MAX_OLIGO_LENGTH,
        MIN_CONST_OLIGO_LENGTH,
    ) -> bool:
        if nd1.index + 3 >= nd2.index:
            return False

        # segments with variable positions between them are required to be at a certain length
        if any([nd1.index < p < nd2.index for p in dna_var_poss]):
            if not (MIN_OLIGO_LENGTH < nd2.index - nd1.index + 3 < MAX_OLIGO_LENGTH):
                return False
        else:
            if nd2.index - nd1.index + 3 < MIN_CONST_OLIGO_LENGTH:
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


# TODO move somewhere
def _find_dna_var_poss(var_poss: List[int]) -> List[int]:
    all_poss: List[int] = []
    for pos in var_poss:
        all_poss.append(pos)
        all_poss.append(pos + 1)
        all_poss.append(pos + 2)
    return all_poss
