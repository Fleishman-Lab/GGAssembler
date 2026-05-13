#!/usr/bin/env python
"""Time the current Python-to-Rust colorful shortest-path wrapper.

Example:
    python benchmarks/bench_colorful_path.py --runs 5 --nodes 300 --edge-probability 0.04
"""

import argparse
import json
import random
import statistics
import time

import networkx as nx
import numpy as np

from dawdlib.dijkstra import colorful, colourful_dijkstra


SEARCH_MODE_FLAGS = {
    "plain": (False, False),
    "astar": (True, False),
    "dominance": (False, True),
    "astar-dominance": (True, True),
}


def build_finder(nodes, edge_probability, colors, seed, resident_graph):
    rng = random.Random(seed)
    graph = nx.DiGraph()
    graph.add_nodes_from(range(nodes))

    for node in range(nodes - 1):
        graph.add_edge(node, node + 1, weight=rng.randint(1, 10))

    for src in range(nodes):
        for dst in range(src + 2, nodes):
            if rng.random() < edge_probability:
                graph.add_edge(src, dst, weight=rng.randint(1, 10))

    finder = object.__new__(colorful.ShortestPathFinder)
    finder.graph = graph
    finder.ggdata = None
    finder.source = 0
    finder.target = nodes - 1
    finder.node_map = {node: node for node in graph.nodes}
    finder.idx_node_map = {node: node for node in graph.nodes}
    finder.graph_edges = [
        (src, dst, data["weight"]) for src, dst, data in graph.edges(data=True)
    ]
    finder.gate_colors = {
        node: [rng.randrange(colors)] if node not in (0, nodes - 1) else []
        for node in graph.nodes
    }
    finder.num_of_colors = [colors]
    finder.all_gate_colors = set(range(colors))
    finder._all_gate_colors_order = list(finder.all_gate_colors)
    finder._gate_color_id_map = {
        color: index for index, color in enumerate(finder._all_gate_colors_order)
    }
    finder._gate_color_ids_by_node = [
        tuple(finder._gate_color_id_map[color] for color in finder.gate_colors[node])
        for node in graph.nodes
    ]
    gate_color_ids_by_node = [
        list(color_ids) for color_ids in finder._gate_color_ids_by_node
    ]
    if resident_graph == "digraph":
        finder._rust_finder = colourful_dijkstra.ColourfulPathFinderDiGraph.with_gate_colors(
            finder.graph_edges,
            finder.node_map[finder.source],
            finder.node_map[finder.target],
            len(finder.node_map),
            gate_color_ids_by_node,
            len(finder._gate_color_id_map),
        )
    else:
        finder._rust_finder = colourful_dijkstra.ColourfulPathFinder.with_gate_colors(
            finder.graph_edges,
            finder.node_map[finder.source],
            finder.node_map[finder.target],
            len(finder.node_map),
            gate_color_ids_by_node,
            len(finder._gate_color_id_map),
        )
    return finder


def time_colab_loop(
    finder, min_gates, max_gates, retries, use_a_star, use_dominance, no_colors
):
    started = time.perf_counter()
    found = 0
    for max_gates_for_run in range(min_gates, max_gates + 1):
        for _ in range(retries):
            path = finder.find_shortest_path(
                len_cutoff=max_gates_for_run,
                no_colors=no_colors or max_gates_for_run + 1,
                use_a_star=use_a_star,
                use_dominance=use_dominance,
            )
            found += bool(path)
    return time.perf_counter() - started, found


def time_find_many(
    finder, min_gates, max_gates, retries, seed, use_a_star, use_dominance, no_colors
):
    started = time.perf_counter()
    results = finder.find_many(
        min_gates,
        max_gates,
        retries,
        seed=seed,
        use_a_star=use_a_star,
        use_dominance=use_dominance,
        no_colors=no_colors,
    )
    return time.perf_counter() - started, len(results)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=int, default=300)
    parser.add_argument("--edge-probability", type=float, default=0.04)
    parser.add_argument("--colors", type=int, default=32)
    parser.add_argument("--min-gates", type=int, default=12)
    parser.add_argument("--max-gates", type=int, default=26)
    parser.add_argument("--retries", type=int, default=20)
    parser.add_argument("--runs", type=int, default=5)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--no-colors",
        type=int,
        default=None,
        help="Override the randomized color bit budget for find-many mode. Default is max_gates + 1 for each gate limit.",
    )
    parser.add_argument(
        "--resident-graph",
        choices=["graphmap", "digraph"],
        default="graphmap",
    )
    parser.add_argument(
        "--mode",
        choices=["legacy-loop", "find-many"],
        default="legacy-loop",
    )
    parser.add_argument(
        "--search-mode",
        choices=sorted(SEARCH_MODE_FLAGS),
        default="astar-dominance",
    )
    args = parser.parse_args()
    use_a_star, use_dominance = SEARCH_MODE_FLAGS[args.search_mode]

    finder = build_finder(
        nodes=args.nodes,
        edge_probability=args.edge_probability,
        colors=args.colors,
        seed=args.seed,
        resident_graph=args.resident_graph,
    )

    samples = []
    found_counts = []
    for run in range(args.runs):
        random.seed(args.seed + run)
        np.random.seed(args.seed + run)
        if args.mode == "find-many":
            elapsed, found = time_find_many(
                finder,
                min_gates=args.min_gates,
                max_gates=args.max_gates,
                retries=args.retries,
                seed=args.seed + run,
                use_a_star=use_a_star,
                use_dominance=use_dominance,
                no_colors=args.no_colors,
            )
        else:
            elapsed, found = time_colab_loop(
                finder,
                min_gates=args.min_gates,
                max_gates=args.max_gates,
                retries=args.retries,
                use_a_star=use_a_star,
                use_dominance=use_dominance,
                no_colors=args.no_colors,
            )
        samples.append(elapsed)
        found_counts.append(found)

    result = {
        "nodes": args.nodes,
        "edges": len(finder.graph_edges),
        "colors": args.colors,
        "min_gates": args.min_gates,
        "max_gates": args.max_gates,
        "retries": args.retries,
        "runs": args.runs,
        "no_colors": args.no_colors,
        "resident_graph": args.resident_graph,
        "mode": args.mode,
        "search_mode": args.search_mode,
        "use_a_star": use_a_star,
        "use_dominance": use_dominance,
        "seconds": samples,
        "mean_seconds": statistics.mean(samples),
        "median_seconds": statistics.median(samples),
        "stdev_seconds": statistics.stdev(samples) if len(samples) > 1 else 0.0,
        "found_counts": found_counts,
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
