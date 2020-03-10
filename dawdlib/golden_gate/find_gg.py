import logging
import os
from itertools import chain
from typing import Callable, Dict, Generator, Iterable, List

import networkx as nx
import pandas as pd

from dawdlib.degenerate_dna.deg_table import TableColNames
from dawdlib.dijkstra.len_limit import all_shortest_paths
from dawdlib.golden_gate.gate import Gate, GateSet
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.gate_restriction import create_path_validator
from dawdlib.golden_gate.graph_maker import GraphMaker, make_default_graph
from dawdlib.golden_gate.utils import Requirements, expand_dna_var_poss, parse_dna

GATESET_FILENAME = "{}_oligos_gate_set_{}.csv"

LOGGER = logging.getLogger(__name__)


def gen_gate_sets(
    graph: nx.Graph,
    source: Gate,
    target: Gate,
    is_valid_path: Callable,
    no_paths: int = 1,
    len_cutoff: int = None,
    weight="weight",
    **kwargs,
) -> Generator[List[Gate], None, None]:
    path_counter = 0
    for gate_set in all_shortest_paths(
        graph, source, target, weight=weight, len_cutoff=len_cutoff, **kwargs
    ):
        if is_valid_path(gate_set):
            yield gate_set
            path_counter += 1
            if path_counter >= no_paths:
                break


def gen_oligo_gateset(
    graph: nx.Graph,
    source: Gate,
    target: Gate,
    is_valid_path: Callable,
    no_paths: int = 1,
    len_cutoff: int = None,
    **kwargs,
) -> Generator[GateSet, None, None]:
    no_oligos = str(len_cutoff) if len_cutoff is not None else "found"
    gate_set_gen = gen_gate_sets(
        graph, source, target, is_valid_path, no_paths, len_cutoff, **kwargs
    )
    for i, gate_set in enumerate(gate_set_gen):
        df = pd.DataFrame(gate_set)
        yield GateSet(no_oligos, i, df)


def gen_gateset_for_oligo_no(
    graph: nx.Graph,
    source: Gate,
    target: Gate,
    is_valid_path: Callable,
    no_paths: int = 1,
    min_oligos: int = None,
    max_oligos: int = None,
    **kwargs,
) -> Iterable[GateSet]:
    if min_oligos is None or max_oligos is None:
        return gen_oligo_gateset(
            graph, source, target, is_valid_path, no_paths, **kwargs
        )
    oligo_gg_sets = (
        gen_oligo_gateset(graph, source, target, is_valid_path, no_paths, i, **kwargs)
        for i in range(min_oligos, max_oligos)
    )
    return chain(*oligo_gg_sets)


def write_oligo_gatesets(gatesets: Iterable[GateSet], outdir: str):
    try:
        for gateset in gatesets:
            gateset_filename = os.path.join(
                outdir, GATESET_FILENAME.format(gateset.no_oligos, gateset.idx)
            )
            gateset.table.to_csv(gateset_filename, index=False)
    except ValueError:
        LOGGER.exception("Failed writing gateset.")
    except nx.NetworkXNoPath:
        LOGGER.exception("Failed writing gateset! No path was found!")


def deg_table_to_dict(deg_table: pd.DataFrame) -> Dict[int, List[str]]:
    keys: List[int] = (deg_table[TableColNames.DNA_POS.value] - 1).tolist()
    values: List[List[str]] = [
        list(filter(len, sub_list))
        for sub_list in deg_table.filter(
            regex=(TableColNames.AMBIGUOUS_CODONS.value + r"[\d]+")
        ).values
    ]
    return dict(zip(keys, values))


def gate_deg_codons(
    dna: str,
    deg_table: pd.DataFrame,
    ggdata: GGData,
    min_var_oligo_length: int,
    max_var_oligo_length: int,
    min_const_oligo_length: int,
    no_solutions: int = 1,
    min_oligos: int = None,
    max_oligos: int = None,
    gate_self_binding_min: int = 2000,
    gate_crosstalk_max: int = 1000,
) -> Iterable[GateSet]:

    is_valid_path = create_path_validator(
        ggdata, gate_self_binding_min, gate_crosstalk_max
    )
    var_poss = deg_table[TableColNames.DNA_POS.value].tolist()
    var_poss = expand_dna_var_poss(var_poss)
    reqs = Requirements(
        min_var_oligo_length,
        max_var_oligo_length,
        min_const_oligo_length,
        gate_self_binding_min,
        gate_crosstalk_max,
    )
    d_graph, src, target = make_default_graph(
        GraphMaker(ggdata), dna, var_poss, deg_table_to_dict(deg_table), reqs
    )
    return gen_gateset_for_oligo_no(
        d_graph, src, target, is_valid_path, no_solutions, min_oligos, max_oligos
    )


def create_goldengates(
    dna_filename: str,
    deg_table_filename: str,
    output_dir: str,
    min_var_oligo_length: int,
    max_var_oligo_length: int,
    min_const_oligo_length: int,
    no_solutions: int = 1,
    min_oligos: int = None,
    max_oligos: int = None,
    gate_self_binding_min: int = 2000,
    gate_crosstalk_max: int = 1000,
):
    dna = parse_dna(dna_filename)
    deg_table = pd.read_csv(
        deg_table_filename, index_col=False, na_filter=True, keep_default_na=False
    )
    ggdata = GGData()
    gatesets = gate_deg_codons(
        dna,
        deg_table,
        ggdata,
        min_var_oligo_length,
        max_var_oligo_length,
        min_const_oligo_length,
        no_solutions,
        min_oligos,
        max_oligos,
        gate_self_binding_min,
        gate_crosstalk_max,
    )
    write_oligo_gatesets(gatesets, output_dir)
