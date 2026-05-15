import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

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
from benchmarks.bench_colorful_real_data import compact_summary


RUN_REAL_DATA_SEARCH = os.environ.get("GGASSEMBLER_RUN_REAL_DATA_SEARCH") == "1"


ROOT = Path(__file__).resolve().parents[3]
EXAMPLE_DIR = ROOT / "example"
FASTA_PATH = EXAMPLE_DIR / "wt_dna.fasta"
RESFILE_PATH = EXAMPLE_DIR / "input.resfile"
DEG_TABLE_PATH = EXAMPLE_DIR / "deg_table.csv"
FIXED_COLORMAP_PATH = (
    Path(__file__).parent
    / "fixtures"
    / "real_data_fixed_colormap_gates_15_seed_1500001.json"
)


def _example_requirements():
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


@pytest.fixture(scope="module")
def real_example_graph():
    resfile = parse_resfile(str(RESFILE_PATH))
    deg_table = pd.read_csv(DEG_TABLE_PATH, na_filter=True, keep_default_na=False)
    assert set(resfile) == set(deg_table[TableColNames.AA_POS.value].tolist())

    reqs = _example_requirements()
    dna = parse_dna(str(FASTA_PATH)).upper()
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


def _path_summary(graph, ggdata, path):
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
        "cost": cost,
        "fidelity": fidelity,
    }


def _assert_valid_path(graph, ggdata, source, target, path, max_gates):
    assert path[0] == source
    assert path[-1] == target
    for n1, n2 in zip(path[:-1], path[1:]):
        assert graph.has_edge(n1, n2)
    summary = _path_summary(graph, ggdata, path)
    assert summary["gate_count"] <= max_gates
    assert summary["cost"] >= 0


def test_real_data_benchmark_summary_handles_no_found_paths():
    assert compact_summary([]) == {
        "found_count": 0,
        "best_cost": None,
        "best_fidelity": None,
        "shortest_gate_count": None,
    }


def test_real_example_graph_builds_from_checked_in_inputs(real_example_graph):
    graph, _ggdata, source, target = real_example_graph

    assert graph.number_of_nodes() == 609
    assert graph.number_of_edges() == 19853
    assert source in graph
    assert target in graph


def test_real_example_dense_node_order_is_topological(real_example_graph):
    graph, ggdata, source, target = real_example_graph
    finder = colorful.ShortestPathFinder(graph, ggdata, source, target)

    assert finder.node_map[source] == 0
    assert finder.node_map[target] == len(finder.node_map) - 1
    assert all(src < dst for src, dst, _weight in finder.graph_edges)


def test_real_example_digraph_matches_fixed_colormap(real_example_graph):
    graph, ggdata, source, target = real_example_graph
    finder = colorful.ShortestPathFinder(graph, ggdata, source, target)
    fixed_colormap = json.loads(FIXED_COLORMAP_PATH.read_text())

    assert fixed_colormap["node_count"] == len(finder.node_map)
    assert fixed_colormap["edge_count"] == len(finder.graph_edges)
    assert fixed_colormap["source_index"] == finder.node_map[source]
    assert fixed_colormap["target_index"] == finder.node_map[target]

    node_colors = list(fixed_colormap["mask_by_node_index"])
    max_gates = fixed_colormap["max_gates"]
    digraph_finder = colourful_dijkstra.ColourfulPathFinderDiGraph(
        finder.graph_edges,
        finder.node_map[source],
        finder.node_map[target],
        len(finder.node_map),
    )

    digraph_path = digraph_finder.find_shortest_path_with_node_colors(
        node_colors, max_gates
    )

    if not digraph_path:
        return

    original_path = [finder.idx_node_map[node_idx] for node_idx in digraph_path]
    _assert_valid_path(graph, ggdata, source, target, original_path, max_gates)


def test_real_example_dag_dp_fixed_colormap_returns_valid_path(real_example_graph):
    graph, ggdata, source, target = real_example_graph
    finder = colorful.ShortestPathFinder(graph, ggdata, source, target)
    fixed_colormap = json.loads(FIXED_COLORMAP_PATH.read_text())
    node_colors = list(fixed_colormap["mask_by_node_index"])
    max_gates = fixed_colormap["max_gates"]

    rust_finder = colourful_dijkstra.ColourfulPathFinderDiGraph(
        finder.graph_edges,
        finder.node_map[source],
        finder.node_map[target],
        len(finder.node_map),
    )

    dag_path = rust_finder.find_shortest_path_dag_dp_with_node_colors(
        node_colors, max_gates
    )

    assert dag_path

    dag_original = [finder.idx_node_map[idx] for idx in dag_path]
    dag_summary = _path_summary(graph, ggdata, dag_original)

    assert dag_summary["cost"] == 2035.0
    _assert_valid_path(graph, ggdata, source, target, dag_original, max_gates)


@pytest.mark.skipif(
    not RUN_REAL_DATA_SEARCH,
    reason="set GGASSEMBLER_RUN_REAL_DATA_SEARCH=1 to run real-data colorful search",
)
def test_real_example_colorful_path_invariants(real_example_graph):
    graph, ggdata, source, target = real_example_graph
    finder = colorful.ShortestPathFinder(graph, ggdata, source, target)
    found_paths = defaultdict(list)
    seed = 1
    min_gates = int(os.environ.get("GGASSEMBLER_REAL_DATA_MIN_GATES", "26"))
    max_gates = int(os.environ.get("GGASSEMBLER_REAL_DATA_MAX_GATES", "26"))
    retries = int(os.environ.get("GGASSEMBLER_REAL_DATA_RETRIES", "1"))

    for max_gates_for_run in range(min_gates, max_gates + 1):
        for retry_index in range(retries):
            np.random.seed(seed + retry_index + max_gates_for_run * 100_000)
            path = finder.find_shortest_path(
                len_cutoff=max_gates_for_run,
                no_colors=max_gates_for_run + 1,
            )
            if not path:
                continue
            _assert_valid_path(
                graph, ggdata, source, target, path, max_gates_for_run
            )
            found_paths[max_gates_for_run].append(
                _path_summary(graph, ggdata, path)
            )

    assert found_paths, "fixed-seed real example run found no colorful paths"
