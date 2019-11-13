from typing import Generator, List, Set

import networkx as nx

from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData


def find_compatible_paths(
    d_graph: nx.DiGraph,
    source: Gate,
    depth_limit: int,
    gg_data: GGData,
    restriction_graph: nx.Graph,
) -> Generator[List[Gate], None, None]:
    nodes = [source]
    visited: Set[Gate] = set()
    for start in nodes:
        if start in visited:
            continue
        visited.add(start)
        stack = [
            (
                start,
                depth_limit,
                iter(d_graph[start]),
                [start],
                set(restriction_graph.neighbors(source)),
            )
        ]
        while stack:
            _, depth_now, children, pth, compatible_neighbors = stack[-1]
            try:
                child = next(children)
                if child not in compatible_neighbors:
                    continue
                if child not in visited:
                    yield pth + [child]
                    visited.add(child)
                    if depth_now > 1:
                        stack.append(
                            (
                                child,
                                depth_now - 1,
                                iter(d_graph[child]),
                                pth + [child],
                                compatible_neighbors
                                ^ set(restriction_graph.neighbors(child)),
                            )
                        )
            except StopIteration:
                stack.pop()


def dijkstra_all_paths(
    d_graph: nx.DiGraph,
    sources: List[Gate],
    sinks: List[Gate],
    depth_limit: int,
    gg_data: GGData,
) -> Generator[List[Gate], None, None]:

    for source in sources:
        for sink in sinks:
            path = nx.dijkstra_path(d_graph, source, sink, "weight")
            if not gg_data.gate_set_has_off_target([a.bps for a in path]):
                yield path
    yield []
