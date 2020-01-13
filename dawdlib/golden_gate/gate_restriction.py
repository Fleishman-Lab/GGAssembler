from itertools import combinations
from typing import Callable, List

import networkx as nx
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gate import Gate


def gen_gate_restriction_graph(
    ggdata: GGData, gate_self_binding_min: int, gate_crosstalk_max: int
) -> nx.Graph:
    nodes = ggdata.filter_self_binding_gates(gate_self_binding_min)
    edges = [
        (v1, v2)
        for v1, v2 in combinations(nodes, 2)
        if ggdata.gates_all_scores(v1, v2) >= gate_crosstalk_max
    ]

    incomp_graph = nx.Graph()
    incomp_graph.add_nodes_from(nodes)
    incomp_graph.add_edges_from(edges)
    return incomp_graph


def create_path_validator(
    ggdata: GGData, gate_self_binding_min: int, gate_crosstalk_max: int
) -> Callable:
    incomp_graph = gen_gate_restriction_graph(
        ggdata, gate_self_binding_min, gate_crosstalk_max
    )

    def valid_path(gg_path: List[Gate]) -> bool:
        return len(incomp_graph.subgraph([node.bps for node in gg_path]).edges) == 0

    return valid_path
