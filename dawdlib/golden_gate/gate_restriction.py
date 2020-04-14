from typing import Callable, List

import networkx as nx

from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData


def gen_gate_restriction_graph(ggdata: GGData) -> nx.Graph:
    nodes = ggdata.filter_self_binding_gates()
    edges = ggdata.restriction_edges(nodes)

    incomp_graph = nx.Graph()
    incomp_graph.add_nodes_from(nodes)
    incomp_graph.add_edges_from(edges)
    incomp_graph.remove_edges_from(nx.selfloop_edges(incomp_graph))

    return incomp_graph


def create_path_validator(ggdata: GGData) -> Callable:
    incomp_graph = gen_gate_restriction_graph(ggdata)

    def valid_path(gg_path: List[Gate]) -> bool:
        return len(incomp_graph.subgraph([node.bps for node in gg_path]).edges) == 0

    return valid_path
