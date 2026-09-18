from itertools import combinations
from typing import Callable, Dict, List, NamedTuple, Optional, Sequence, Tuple, Union

import networkx as nx
import numpy as np

from dawdlib.golden_gate.constants import (
    CONST_COST,
    SINK_NODE,
    SOURCE_NODE,
    VAR_ADD_COST,
)
from dawdlib.golden_gate.gate import Gate, PseudoGate
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.utils import Requirements, syn_muts

SOURCE = PseudoGate(idx=-1, bps=SOURCE_NODE, src_or_target=True)
TARGET = PseudoGate(idx=1, bps=SINK_NODE, src_or_target=True)

INDIVIDUAL_OLIGO_TIERS = ((13, 60, 13), (61, 80, 17))


class ExactPhysicalPart(NamedTuple):
    payload_sequence: str
    ordered_sequence: str
    source_type: str
    is_wt: bool
    cost_cents: int
    primer_length: int
    primer_count: int
    variant_ids: Tuple[str, ...]


class ExactIntervalPlan(NamedTuple):
    parts: Tuple[ExactPhysicalPart, ...]
    variant_part_indices: Tuple[Tuple[str, int], ...]
    cost_cents: int


class GraphMaker:
    def __init__(self, ggdata):
        self.ggdata: GGData = ggdata

    def create_default_valid_edge_func(
        self,
        dna_var_poss: List[int],
        min_oligo_length: int,
        max_oligo_length: int,
        min_const_oligo_length: int,
        min_fidelity: float,
    ) -> Callable[[Gate, Gate], bool]:
        dna_var_arr = np.array(dna_var_poss)

        def _is_valid_edge(nd1: Gate, nd2: Gate) -> bool:
            if nd1.overlap(nd2):
                return False
            # segments with variable positions between them
            # are required to be at a certain length
            start = min(nd1.idx, nd2.idx)
            end = max(nd1.idx, nd2.idx)
            psuedogates = isinstance(nd1, PseudoGate) or isinstance(nd2, PseudoGate)
            if np.any(np.logical_and(start < dna_var_arr, dna_var_arr < end)):
                if not (min_oligo_length < nd2 - nd1 < max_oligo_length):
                    return False
            else:
                if not psuedogates and (nd2 - nd1 < min_const_oligo_length):
                    return False
            if (
                not psuedogates
                and self.ggdata.overhangs_fidelity(nd1.bps, nd2.bps) < min_fidelity
            ):
                return False
            return True

        return _is_valid_edge


def make_nodes(
    d_graph,
    dna: str,
    is_valid_node: Callable[[Gate], bool],
    add_syn: bool = False,
    glength: int = 4,
) -> None:
    if add_syn and len(dna) % 3 == 0:
        add_syn = True
    possible_nodes: List[Gate] = []
    for ind in range(len(dna) - glength - 1):
        gate_idxs = slice(ind, ind + glength)
        fcw = dna[gate_idxs]
        possible_nodes.append(Gate(ind, fcw, False, gatelength=glength))
    if add_syn:
        possible_nodes += syn_muts(dna)
    d_graph.add_nodes_from(filter(is_valid_node, possible_nodes))


def create_default_valid_node_function(
    acceptable_fcws: List[str], var_dna_poss: List[int]
) -> Callable[[Gate], bool]:
    var_dna_arr = np.array(var_dna_poss)

    def _is_valid_node(gate: Gate) -> bool:

        if gate.bps not in acceptable_fcws:
            return False
        gate_span = gate.span()
        if np.any(
            np.logical_and(gate_span[0] <= var_dna_arr, var_dna_arr <= gate_span[1])
        ):
            return False
        return True

    return _is_valid_node


def create_default_weight_func(
    dna_pos_n_codons: Dict[int, List[str]],
    oligo_addition: int = VAR_ADD_COST,
    const_cost: int = CONST_COST,
) -> Callable[[Gate, Gate], int]:
    var_pos_arr = np.array(list(dna_pos_n_codons.keys()))

    def edge_weight(nd1: Gate, nd2: Gate) -> int:
        start = min(nd1.idx, nd2.idx)
        end = max(nd1.idx, nd2.idx)
        psuedogates = isinstance(nd1, PseudoGate) or isinstance(nd2, PseudoGate)
        var_pos = var_pos_arr[np.logical_and(start < var_pos_arr, var_pos_arr < end)]
        if var_pos.size:
            return (nd2 - nd1 + oligo_addition) * np.prod(
                list(len(dna_pos_n_codons[pos]) for pos in var_pos)
            )
        return const_cost if not psuedogates else 0

    return edge_weight


def _gate_oligo_slice(g1: Gate, g2: Gate) -> slice:
    g1_idx = g1.idx if not g1.src_or_target else None
    g2_idx = g2.span()[1] + 1 if not g2.src_or_target else None
    return slice(g1_idx, g2_idx)


def make_length_tier_cost(
    tiers: Sequence[Tuple[int, int, int]],
    fallback: Optional[Callable[[int], Optional[int]]] = None,
) -> Callable[[int], Optional[int]]:
    def _cost(length: int) -> Optional[int]:
        matching_costs = [
            length * price_per_nt
            for min_len, max_len, price_per_nt in tiers
            if min_len <= length <= max_len
        ]
        if matching_costs:
            return min(matching_costs)
        if fallback is not None:
            return fallback(length)
        return None

    return _cost


def _validated_length_cost(
    length_cost: Callable[[int], Optional[int]], length: int
) -> Optional[int]:
    cost = length_cost(length)
    if cost is None:
        return None
    if isinstance(cost, bool) or not isinstance(cost, int) or cost < 0:
        raise ValueError("Oligo length cost must be a non-negative integer or None")
    return cost


def plan_exact_interval(
    nd1: Gate,
    nd2: Gate,
    wt_dna: str,
    variant_dnas: Dict[str, str],
    prefix: str,
    suffix: str,
    segment_length_cost: Callable[[int], Optional[int]],
    min_const_oligo_length: int = 0,
    primer_length: int = 20,
    primer_price_per_nt_cents: int = 13,
    allow_wt_direct_synthesis: bool = False,
) -> Optional[ExactIntervalPlan]:
    o_slice = _gate_oligo_slice(nd1, nd2)
    wt_segment = wt_dna[o_slice]
    terminal = isinstance(nd1, PseudoGate) or isinstance(nd2, PseudoGate)
    ordered_prefix = prefix if not nd1.src_or_target else ""
    ordered_suffix = suffix if not nd2.src_or_target else ""
    variants_by_segment: Dict[str, List[str]] = {}
    for variant_id, variant_dna in variant_dnas.items():
        variants_by_segment.setdefault(variant_dna[o_slice], []).append(variant_id)

    parts = []
    variant_part_indices = []
    pcr_cost = 2 * primer_length * primer_price_per_nt_cents
    for segment, variant_ids in variants_by_segment.items():
        is_wt = segment == wt_segment
        ordered_sequence = ordered_prefix + segment + ordered_suffix
        primer_count = 0
        part_primer_length = 0

        if terminal and is_wt:
            source_type = "terminal_free"
            ordered_sequence = ""
            cost_cents = 0
        elif is_wt:
            synthesis_cost = None
            if allow_wt_direct_synthesis:
                synthesis_cost = _validated_length_cost(
                    segment_length_cost, len(ordered_sequence)
                )
            pcr_is_feasible = len(segment) >= min_const_oligo_length
            if synthesis_cost is not None and (
                not pcr_is_feasible or synthesis_cost < pcr_cost
            ):
                source_type = "individual_oligo"
                cost_cents = synthesis_cost
            elif pcr_is_feasible:
                source_type = "wt_pcr"
                ordered_sequence = ""
                cost_cents = pcr_cost
                primer_count = 2
                part_primer_length = primer_length
            else:
                return None
        else:
            cost_cents = _validated_length_cost(
                segment_length_cost, len(ordered_sequence)
            )
            if cost_cents is None:
                return None
            source_type = "individual_oligo"

        part_index = len(parts)
        parts.append(
            ExactPhysicalPart(
                payload_sequence=segment,
                ordered_sequence=ordered_sequence,
                source_type=source_type,
                is_wt=is_wt,
                cost_cents=cost_cents,
                primer_length=part_primer_length,
                primer_count=primer_count,
                variant_ids=tuple(variant_ids),
            )
        )
        variant_part_indices.extend(
            (variant_id, part_index) for variant_id in variant_ids
        )

    return ExactIntervalPlan(
        parts=tuple(parts),
        variant_part_indices=tuple(variant_part_indices),
        cost_cents=sum(part.cost_cents for part in parts),
    )


def create_exact_variant_weight_func(
    wt_dna: str,
    variant_dnas: Dict[str, str],
    prefix: str = "",
    suffix: str = "",
    segment_length_cost: Optional[Callable[[int], Optional[int]]] = None,
    min_const_oligo_length: int = 0,
    primer_length: int = 20,
    primer_price_per_nt_cents: int = 13,
    allow_wt_direct_synthesis: bool = False,
) -> Callable[[Gate, Gate], Optional[int]]:
    if segment_length_cost is None:
        segment_length_cost = make_length_tier_cost(INDIVIDUAL_OLIGO_TIERS)

    def edge_weight(nd1: Gate, nd2: Gate) -> Optional[int]:
        plan = plan_exact_interval(
            nd1,
            nd2,
            wt_dna,
            variant_dnas,
            prefix,
            suffix,
            segment_length_cost,
            min_const_oligo_length,
            primer_length,
            primer_price_per_nt_cents,
            allow_wt_direct_synthesis,
        )
        return None if plan is None else plan.cost_cents

    return edge_weight


def make_edges(
    d_graph: nx.Graph,
    is_valid_edge: Callable[[Gate, Gate], bool],
    edge_weight: Callable[[Gate, Gate], Union[float, int]],
) -> None:
    d_graph.add_weighted_edges_from(
        map(
            lambda x: x + (edge_weight(x[0], x[1]),),
            filter(
                lambda x: x[0].idx < x[1].idx and is_valid_edge(x[0], x[1]),
                combinations(sorted(d_graph.nodes, key=lambda x: x.idx), 2),
            ),
        )
    )


def build_custom_graph(
    dna: str,
    is_valid_node: Callable[[Gate], bool],
    is_valid_edge: Callable[[Gate, Gate], bool],
    edge_weight: Callable[[Gate, Gate], Union[float, int]],
    gate_length: int = 4,
) -> Tuple[nx.Graph, Gate, Gate]:
    d_graph: nx.DiGraph = nx.DiGraph()
    src = SOURCE
    target = TARGET._replace(idx=TARGET.idx + len(dna))

    make_nodes(d_graph, dna, is_valid_node, glength=gate_length)
    d_graph.add_node(src)
    d_graph.add_node(target)
    make_edges(d_graph, is_valid_edge, edge_weight)
    return d_graph, src, target


def make_default_graph(
    gm: GraphMaker,
    dna: str,
    var_poss: List[int],
    dna_pos_n_codons: Dict[int, List[str]],
    reqs: Requirements,
    gatelength: int = 4,
) -> Tuple[nx.Graph, Gate, Gate]:
    acceptable_fcws = gm.ggdata.filter_self_binding_gates(reqs.filter_gc_overhangs)
    is_valid_node = create_default_valid_node_function(acceptable_fcws, var_poss)
    is_valid_edge = gm.create_default_valid_edge_func(
        var_poss,
        reqs.min_oligo_length,
        reqs.max_oligo_length,
        reqs.min_const_oligo_length,
        reqs.min_fidelity,
    )
    edge_weight = create_default_weight_func(
        dna_pos_n_codons, reqs.oligo_addition, reqs.const_cost
    )
    return build_custom_graph(
        dna, is_valid_node, is_valid_edge, edge_weight, gate_length=gatelength
    )


def make_no_degenerate_graph(
    gm: GraphMaker,
    dna: str,
    reqs: Requirements,
    gatelength: int = 4,
) -> Tuple[nx.Graph, Gate, Gate]:
    return make_default_graph(gm, dna, [], {}, reqs, gatelength)


def make_exact_variant_graph(
    gm: GraphMaker,
    dna: str,
    variant_dnas: Dict[str, str],
    reqs: Requirements,
    segment_length_cost: Optional[Callable[[int], Optional[int]]] = None,
    primer_length: int = 20,
    primer_price_per_nt_cents: int = 13,
    allow_wt_direct_synthesis: bool = False,
    overhang_filter: Optional[Callable[[str], bool]] = None,
    gatelength: int = 4,
) -> Tuple[nx.Graph, Gate, Gate]:
    if not variant_dnas:
        raise ValueError("At least one exact variant is required")
    dna = dna.upper()
    variant_dnas = {
        variant_id: variant_dna.upper()
        for variant_id, variant_dna in variant_dnas.items()
    }
    if any(len(variant_dna) != len(dna) for variant_dna in variant_dnas.values()):
        raise ValueError("Every exact variant DNA must have the WT DNA length")
    if any(
        set(variant_dna.upper()) - set("ACGT") for variant_dna in variant_dnas.values()
    ):
        raise ValueError("Exact variant DNA must contain only A, C, G, and T")
    if set(dna.upper()) - set("ACGT"):
        raise ValueError("WT DNA must contain only A, C, G, and T")

    if segment_length_cost is None:
        segment_length_cost = make_length_tier_cost(INDIVIDUAL_OLIGO_TIERS)
    acceptable_fcws = set(gm.ggdata.filter_self_binding_gates(reqs.filter_gc_overhangs))

    def is_valid_node(gate: Gate) -> bool:
        return (
            gate.bps in acceptable_fcws
            and (overhang_filter is None or overhang_filter(gate.bps))
            and all(
                variant_dna[gate.idx : gate.idx + gate.gatelength] == gate.bps
                for variant_dna in variant_dnas.values()
            )
        )

    plans: Dict[Tuple[Gate, Gate], ExactIntervalPlan] = {}

    def is_valid_edge(nd1: Gate, nd2: Gate) -> bool:
        if nd1.overlap(nd2):
            return False
        pseudogates = isinstance(nd1, PseudoGate) or isinstance(nd2, PseudoGate)
        if (
            not pseudogates
            and gm.ggdata.overhangs_fidelity(nd1.bps, nd2.bps) < reqs.min_fidelity
        ):
            return False
        plan = plan_exact_interval(
            nd1,
            nd2,
            dna,
            variant_dnas,
            reqs.oligo_prefix,
            reqs.oligo_suffix,
            segment_length_cost,
            reqs.min_const_oligo_length,
            primer_length,
            primer_price_per_nt_cents,
            allow_wt_direct_synthesis,
        )
        if plan is None:
            return False
        plans[(nd1, nd2)] = plan
        return True

    def edge_weight(nd1: Gate, nd2: Gate) -> int:
        return plans[(nd1, nd2)].cost_cents

    graph, source, target = build_custom_graph(
        dna, is_valid_node, is_valid_edge, edge_weight, gate_length=gatelength
    )
    nx.set_edge_attributes(
        graph,
        {(nd1, nd2): plan for (nd1, nd2), plan in plans.items()},
        "exact_plan",
    )
    return graph, source, target
