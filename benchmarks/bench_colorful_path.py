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

from dawdlib.dijkstra import colorful


def build_finder(nodes, edge_probability, colors, seed):
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
    return finder


def time_colab_loop(finder, min_gates, max_gates, retries):
    started = time.perf_counter()
    found = 0
    for max_gates_for_run in range(min_gates, max_gates + 1):
        for _ in range(retries):
            path = finder.find_shortest_path(
                len_cutoff=max_gates_for_run,
                no_colors=max_gates_for_run + 1,
            )
            found += bool(path)
    return time.perf_counter() - started, found


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
    args = parser.parse_args()

    finder = build_finder(
        nodes=args.nodes,
        edge_probability=args.edge_probability,
        colors=args.colors,
        seed=args.seed,
    )

    samples = []
    found_counts = []
    for run in range(args.runs):
        random.seed(args.seed + run)
        np.random.seed(args.seed + run)
        elapsed, found = time_colab_loop(
            finder,
            min_gates=args.min_gates,
            max_gates=args.max_gates,
            retries=args.retries,
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
        "seconds": samples,
        "mean_seconds": statistics.mean(samples),
        "median_seconds": statistics.median(samples),
        "stdev_seconds": statistics.stdev(samples) if len(samples) > 1 else 0.0,
        "found_counts": found_counts,
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
