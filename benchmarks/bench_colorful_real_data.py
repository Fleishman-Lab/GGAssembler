#!/usr/bin/env python
"""Time colorful shortest path on the checked-in Golden Gate example data.

The benchmark intentionally loads the checked-in degenerate table instead of
calling generate_deg_csv, because generate_deg_csv may fetch codon usage data
from the network through synbiochem.
"""

import argparse
import json
import statistics
import time
from pathlib import Path

import numpy as np
import pandas as pd

from dawdlib.degenerate_dna.deg_table import TableColNames
from dawdlib.dijkstra import colorful, colourful_dijkstra
from dawdlib.golden_gate.find_gg import deg_table_to_dict
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.graph_maker import GraphMaker, make_default_graph
from dawdlib.golden_gate.utils import (
    RequirementsFactory,
    expand_dna_var_poss,
    parse_dna,
    parse_resfile,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FASTA = ROOT / "example" / "wt_dna.fasta"
DEFAULT_RESFILE = ROOT / "example" / "input.resfile"
DEFAULT_DEG_TABLE = ROOT / "example" / "deg_table.csv"
DEFAULT_FIXED_COLORMAP = (
    ROOT
    / "dawdlib"
    / "dijkstra"
    / "tests"
    / "fixtures"
    / "real_data_fixed_colormap_gates_15_seed_1500001.json"
)


def example_requirements():
    return RequirementsFactory(
        min_oligo_length=4,
        max_oligo_length=100,
        min_const_oligo_length=15,
        min_efficiency=0.25,
        min_fidelity=0.1,
        oligo_prefix="GACATTGGTCTCA",
        oligo_suffix="TGAGACCAACGACGCCGTACTCTTTGTCAAC",
        const_cost=40,
        filter_gc_overhangs=False,
    )


def build_real_example(fasta_path, resfile_path, deg_table_path):
    resfile = parse_resfile(str(resfile_path))
    deg_table = pd.read_csv(deg_table_path, na_filter=True, keep_default_na=False)
    deg_positions = set(deg_table[TableColNames.AA_POS.value].tolist())
    if set(resfile) != deg_positions:
        raise ValueError("resfile positions do not match degenerate table positions")

    reqs = example_requirements()
    dna = parse_dna(str(fasta_path)).upper()
    ggdata = GGData(
        temperature=reqs.gg_temp,
        hours=reqs.gg_hours,
        min_efficiency=reqs.min_efficiency,
        min_fidelity=reqs.min_fidelity,
    )
    var_poss = expand_dna_var_poss(deg_table[TableColNames.DNA_POS.value].tolist())
    graph, source, target = make_default_graph(
        GraphMaker(ggdata),
        dna,
        var_poss,
        deg_table_to_dict(deg_table),
        reqs,
        4,
    )
    return graph, ggdata, source, target


def summarize_path(graph, ggdata, path):
    inner_path = [node for node in path if not node.src_or_target]
    overhangs = [node.bps for node in inner_path]
    cost = sum(
        graph.edges[n1, n2]["weight"] for n1, n2 in zip(path[:-1], path[1:])
    )
    try:
        fidelity = ggdata.reaction_fidelity(*overhangs)[0]
    except ValueError:
        fidelity = None
    return {
        "gate_count": len(inner_path),
        "cost": float(cost),
        "fidelity": float(fidelity) if fidelity is not None else None,
    }


def summarize_dense_path(finder, graph, ggdata, dense_path):
    if not dense_path:
        return {
            "found": False,
            "dense_path": [],
            "dense_path_length": 0,
            "cost": None,
            "fidelity": None,
            "gate_count": None,
        }
    path = [finder.idx_node_map[node_idx] for node_idx in dense_path]
    summary = summarize_path(graph, ggdata, path)
    return {
        "found": True,
        "dense_path": dense_path,
        "dense_path_length": len(dense_path),
        **summary,
    }


def load_fixed_colormap(path, finder):
    fixed_colormap = json.loads(path.read_text())
    expected = {
        "node_count": len(finder.node_map),
        "edge_count": len(finder.graph_edges),
        "source_index": finder.node_map[finder.source],
        "target_index": finder.node_map[finder.target],
    }
    for key, value in expected.items():
        if fixed_colormap[key] != value:
            raise ValueError(
                f"fixed colormap {key}={fixed_colormap[key]} does not match {value}"
            )
    if len(fixed_colormap["mask_by_node_index"]) != len(finder.node_map):
        raise ValueError("fixed colormap mask count does not match node count")
    return fixed_colormap


def node_color_list_from_fixed_colormap(fixed_colormap):
    return list(fixed_colormap["mask_by_node_index"])


def build_resident_finder(finder, resident_graph):
    if resident_graph == "digraph":
        return colourful_dijkstra.ColourfulPathFinderDiGraph(
            finder.graph_edges,
            finder.node_map[finder.source],
            finder.node_map[finder.target],
            len(finder.node_map),
        )
    return colourful_dijkstra.ColourfulPathFinder(
        finder.graph_edges,
        finder.node_map[finder.source],
        finder.node_map[finder.target],
    )


def build_resident_finder_with_gate_colors(finder, resident_graph):
    gate_color_ids_by_node = [
        list(color_ids) for color_ids in finder._gate_color_ids_by_node
    ]
    if resident_graph == "digraph":
        return colourful_dijkstra.ColourfulPathFinderDiGraph.with_gate_colors(
            finder.graph_edges,
            finder.node_map[finder.source],
            finder.node_map[finder.target],
            len(finder.node_map),
            gate_color_ids_by_node,
            len(finder._gate_color_id_map),
        )
    return colourful_dijkstra.ColourfulPathFinder.with_gate_colors(
        finder.graph_edges,
        finder.node_map[finder.source],
        finder.node_map[finder.target],
        len(finder.node_map),
        gate_color_ids_by_node,
        len(finder._gate_color_id_map),
    )


def time_legacy_loop(finder, graph, ggdata, min_gates, max_gates, retries, seed, run):
    started = time.perf_counter()
    summaries = []
    for max_gates_for_run in range(min_gates, max_gates + 1):
        for retry_index in range(retries):
            np.random.seed(
                seed + run * 1_000_000 + max_gates_for_run * 100_000 + retry_index
            )
            path = finder.find_shortest_path(
                len_cutoff=max_gates_for_run,
                no_colors=max_gates_for_run + 1,
            )
            if path:
                summaries.append(summarize_path(graph, ggdata, path))
    return time.perf_counter() - started, summaries


def time_fixed_colormap(finder, graph, ggdata, fixed_colormap, resident_graph):
    resident_finder = build_resident_finder(finder, resident_graph)
    node_colors = node_color_list_from_fixed_colormap(fixed_colormap)
    max_gates = fixed_colormap["max_gates"]

    started = time.perf_counter()
    dense_path = resident_finder.find_shortest_path_with_node_colors(
        node_colors, max_gates
    )
    elapsed = time.perf_counter() - started
    return elapsed, summarize_dense_path(finder, graph, ggdata, dense_path)


def time_fixed_colormap_comparison(finder, graph, ggdata, fixed_colormap):
    node_colors = node_color_list_from_fixed_colormap(fixed_colormap)
    max_gates = fixed_colormap["max_gates"]
    results = {}
    for resident_graph in ("graphmap", "digraph"):
        resident_finder = build_resident_finder(finder, resident_graph)
        started = time.perf_counter()
        dense_path = resident_finder.find_shortest_path_with_node_colors(
            node_colors, max_gates
        )
        elapsed = time.perf_counter() - started
        results[resident_graph] = {
            "seconds": elapsed,
            **summarize_dense_path(finder, graph, ggdata, dense_path),
        }

    graphmap = results["graphmap"]
    digraph = results["digraph"]
    comparison = {
        "both_found_equal": graphmap["found"] == digraph["found"],
        "path_lengths_equal": graphmap["dense_path_length"]
        == digraph["dense_path_length"],
        "costs_equal": graphmap["cost"] == digraph["cost"],
        "same_dense_path": graphmap["dense_path"] == digraph["dense_path"],
    }
    return results, comparison


def time_find_many(finder, graph, ggdata, min_gates, max_gates, retries, seed):
    if not hasattr(finder, "find_many"):
        raise RuntimeError("ShortestPathFinder.find_many is not implemented yet")

    started = time.perf_counter()
    raw_results = finder.find_many(min_gates, max_gates, retries, seed=seed)
    summaries = [
        summarize_path(graph, ggdata, path)
        for _max_gates, _retry_index, path in raw_results
    ]
    return time.perf_counter() - started, summaries


def compact_summary(all_summaries):
    valid_fidelities = [
        summary["fidelity"]
        for summary in all_summaries
        if summary["fidelity"] is not None
    ]
    best_cost = min((summary["cost"] for summary in all_summaries), default=None)
    return {
        "found_count": len(all_summaries),
        "best_cost": float(best_cost) if best_cost is not None else None,
        "best_fidelity": float(max(valid_fidelities)) if valid_fidelities else None,
        "shortest_gate_count": min(
            (summary["gate_count"] for summary in all_summaries), default=None
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--resfile", type=Path, default=DEFAULT_RESFILE)
    parser.add_argument("--deg-table", type=Path, default=DEFAULT_DEG_TABLE)
    parser.add_argument("--min-gates", type=int, default=14)
    parser.add_argument("--max-gates", type=int, default=15)
    parser.add_argument("--retries", type=int, default=1)
    parser.add_argument("--runs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--fixed-colormap", type=Path, default=DEFAULT_FIXED_COLORMAP)
    parser.add_argument(
        "--resident-graph",
        choices=["graphmap", "digraph"],
        default="graphmap",
    )
    parser.add_argument(
        "--mode",
        choices=[
            "legacy-loop",
            "find-many",
            "fixed-colormap",
            "compare-fixed-colormap",
        ],
        default="legacy-loop",
    )
    args = parser.parse_args()

    setup_started = time.perf_counter()
    graph, ggdata, source, target = build_real_example(
        args.fasta, args.resfile, args.deg_table
    )
    finder = colorful.ShortestPathFinder(graph, ggdata, source, target)
    if args.resident_graph == "digraph":
        finder._rust_finder = build_resident_finder_with_gate_colors(
            finder, args.resident_graph
        )
    setup_seconds = time.perf_counter() - setup_started

    if args.mode in {"fixed-colormap", "compare-fixed-colormap"}:
        fixed_colormap = load_fixed_colormap(args.fixed_colormap, finder)
        if args.mode == "compare-fixed-colormap":
            comparison_samples = []
            for _run in range(args.runs):
                comparison_samples.append(
                    time_fixed_colormap_comparison(
                        finder, graph, ggdata, fixed_colormap
                    )
                )
            graphmap_seconds = [
                sample[0]["graphmap"]["seconds"] for sample in comparison_samples
            ]
            digraph_seconds = [
                sample[0]["digraph"]["seconds"] for sample in comparison_samples
            ]
            latest_results, latest_comparison = comparison_samples[-1]
            result = {
                "input_fasta": str(args.fasta),
                "input_resfile": str(args.resfile),
                "input_deg_table": str(args.deg_table),
                "fixed_colormap": str(args.fixed_colormap),
                "fixed_colormap_seed": fixed_colormap["seed"],
                "nodes": graph.number_of_nodes(),
                "edges": graph.number_of_edges(),
                "mode": args.mode,
                "runs": args.runs,
                "setup_seconds": setup_seconds,
                "graphmap_seconds": graphmap_seconds,
                "graphmap_mean_seconds": statistics.mean(graphmap_seconds),
                "graphmap_median_seconds": statistics.median(graphmap_seconds),
                "digraph_seconds": digraph_seconds,
                "digraph_mean_seconds": statistics.mean(digraph_seconds),
                "digraph_median_seconds": statistics.median(digraph_seconds),
                "latest_results": latest_results,
                "latest_comparison": latest_comparison,
            }
            print(json.dumps(result, indent=2, sort_keys=True))
            return

        samples = []
        summaries = []
        for _run in range(args.runs):
            elapsed, summary = time_fixed_colormap(
                finder, graph, ggdata, fixed_colormap, args.resident_graph
            )
            samples.append(elapsed)
            summaries.append(summary)
        result = {
            "input_fasta": str(args.fasta),
            "input_resfile": str(args.resfile),
            "input_deg_table": str(args.deg_table),
            "fixed_colormap": str(args.fixed_colormap),
            "fixed_colormap_seed": fixed_colormap["seed"],
            "nodes": graph.number_of_nodes(),
            "edges": graph.number_of_edges(),
            "max_gates": fixed_colormap["max_gates"],
            "no_colors": fixed_colormap["no_colors"],
            "mode": args.mode,
            "resident_graph": args.resident_graph,
            "runs": args.runs,
            "setup_seconds": setup_seconds,
            "seconds": samples,
            "mean_seconds": statistics.mean(samples),
            "median_seconds": statistics.median(samples),
            "stdev_seconds": statistics.stdev(samples) if len(samples) > 1 else 0.0,
            "latest_summary": summaries[-1],
        }
        print(json.dumps(result, indent=2, sort_keys=True))
        return

    samples = []
    all_run_summaries = []
    for run in range(args.runs):
        if args.mode == "legacy-loop":
            elapsed, summaries = time_legacy_loop(
                finder,
                graph,
                ggdata,
                args.min_gates,
                args.max_gates,
                args.retries,
                args.seed,
                run,
            )
        else:
            elapsed, summaries = time_find_many(
                finder,
                graph,
                ggdata,
                args.min_gates,
                args.max_gates,
                args.retries,
                args.seed + run,
            )
        samples.append(elapsed)
        all_run_summaries.extend(summaries)

    result = {
        "input_fasta": str(args.fasta),
        "input_resfile": str(args.resfile),
        "input_deg_table": str(args.deg_table),
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "min_gates": args.min_gates,
        "max_gates": args.max_gates,
        "retries": args.retries,
        "runs": args.runs,
        "seed": args.seed,
        "mode": args.mode,
        "resident_graph": args.resident_graph,
        "setup_seconds": setup_seconds,
        "seconds": samples,
        "mean_seconds": statistics.mean(samples),
        "median_seconds": statistics.median(samples),
        "stdev_seconds": statistics.stdev(samples) if len(samples) > 1 else 0.0,
    }
    result.update(compact_summary(all_run_summaries))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
