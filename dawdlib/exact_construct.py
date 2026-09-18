import json
from typing import Any, Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple

import pandas as pd
from Bio.Seq import Seq

from dawdlib.degenerate_dna.codon_selector import CodonSelector
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.utils import check_for_restriction_sites, find_dna_var_poss


class Mutation(NamedTuple):
    aa_pos: int
    from_aa: str
    to_aa: str


class Variant(NamedTuple):
    variant_id: str
    target_aas: str
    mutations: Sequence[Mutation]
    total_score: Optional[float] = None


class ExactConstructOutput(NamedTuple):
    parts: pd.DataFrame
    bom: pd.DataFrame
    summary: Dict[str, Any]


def read_sequence_space_resfile(path: str) -> List[Tuple[int, str]]:
    design_space = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line in {"nataa", "start"}:
                continue
            parts = line.split()
            if len(parts) != 4 or parts[2] != "PIKAA":
                raise ValueError(f"Expected '<pos> <chain> PIKAA <aas>': {line}")
            design_space.append((int(parts[0]), parts[3]))
    return design_space


def read_sequence_space_variants(
    variants_csv_path: str,
    resfile_path: str,
    wt_dna: str,
) -> List[Variant]:
    design_space = read_sequence_space_resfile(resfile_path)
    positions = [pos for pos, _ in design_space]
    allowed_aas = [aas for _, aas in design_space]
    wt_aas = str(Seq(wt_dna).translate())
    df = pd.read_csv(variants_csv_path)

    variants = []
    seen_ids = set()
    for row in df.itertuples(index=False):
        target_aas = str(row.sequence_name).strip()
        if target_aas in seen_ids:
            raise ValueError(f"Duplicate sequence_name: {target_aas}")
        seen_ids.add(target_aas)
        if len(target_aas) != len(positions):
            raise ValueError(
                f"{target_aas} has {len(target_aas)} amino acids, "
                f"expected {len(positions)}"
            )

        mutations = []
        for aa, aa_pos, allowed in zip(target_aas, positions, allowed_aas):
            if aa not in allowed:
                raise ValueError(f"{aa} is not allowed at position {aa_pos}")
            if aa_pos > len(wt_aas):
                raise ValueError(
                    f"Position {aa_pos} exceeds translated WT length {len(wt_aas)}"
                )
            wt_aa = wt_aas[aa_pos - 1]
            if wt_aa not in allowed:
                raise ValueError(f"WT {wt_aa} is not allowed at position {aa_pos}")
            if aa != wt_aa:
                mutations.append(Mutation(aa_pos, wt_aa, aa))

        total_score = getattr(row, "total_score", None)
        if pd.isna(total_score):
            total_score = None
        variants.append(
            Variant(
                variant_id=target_aas,
                target_aas=target_aas,
                mutations=tuple(mutations),
                total_score=total_score,
            )
        )

    return variants


def variant_mutation_positions(variants: Iterable[Variant]) -> List[int]:
    aa_positions = sorted(
        {mutation.aa_pos for variant in variants for mutation in variant.mutations}
    )
    return find_dna_var_poss(aa_positions)


def choose_exact_codon(codon_selector: CodonSelector, amino_acid: str) -> str:
    selected = codon_selector.optimise_codons([amino_acid])
    codons = [
        codon
        for selected_aa in selected.amino_acids
        if selected_aa.amino_acid == amino_acid
        for codon in selected_aa.codons
    ]
    if codons:
        return max(codons, key=lambda codon: (codon.cai, codon.probability)).codon

    for codon in selected.ambiguous_codons:
        if all(base in "ACGT" for base in codon):
            return codon
    raise ValueError(f"Expected exact codon for {amino_acid}, got {selected}")


def build_exact_variant_dnas(
    wt_dna: str,
    variants: Sequence[Variant],
    codon_selector: CodonSelector,
    enzymes: Optional[Sequence[str]] = None,
) -> Dict[str, str]:
    wt_dna = wt_dna.upper()
    variant_ids = [variant.variant_id for variant in variants]
    if len(variant_ids) != len(set(variant_ids)):
        raise ValueError("Exact variant IDs must be unique")
    if set(wt_dna) - set("ACGT"):
        raise ValueError("WT DNA must contain only A, C, G, and T")
    result = {}

    for variant in variants:
        dna = list(wt_dna)
        for mutation in variant.mutations:
            dna_idx = (mutation.aa_pos - 1) * 3
            codon = choose_exact_codon(codon_selector, mutation.to_aa)
            dna[dna_idx : dna_idx + 3] = codon

        variant_dna = "".join(dna)
        translated = str(Seq(variant_dna).translate())
        for mutation in variant.mutations:
            if translated[mutation.aa_pos - 1] != mutation.to_aa:
                raise ValueError(
                    f"{variant.variant_id} does not encode {mutation.to_aa} "
                    f"at position {mutation.aa_pos}"
                )
        if enzymes:
            ok, enzyme, sites = check_for_restriction_sites(variant_dna, list(enzymes))
            if not ok:
                raise ValueError(
                    f"{variant.variant_id} introduces {enzyme} sites at {sites}"
                )
        result[variant.variant_id] = variant_dna

    return result


def exact_variant_tables(
    path: Sequence[Gate],
    graph,
    variant_dnas: Dict[str, str],
    fidelity_metrics: Optional[Dict[str, float]] = None,
) -> ExactConstructOutput:
    if not variant_dnas:
        raise ValueError("At least one exact variant is required")

    part_rows = []
    bom_rows = []
    for interval_number, (g1, g2) in enumerate(zip(path[:-1], path[1:]), 1):
        plan = graph.edges[g1, g2].get("exact_plan")
        if plan is None:
            raise ValueError(f"Missing exact interval plan for {g1.idx}-{g2.idx}")

        interval_id = f"I{interval_number}"
        gate1 = json.dumps(g1)
        gate2 = json.dumps(g2)
        part_ids = []
        for part_number, part in enumerate(plan.parts, 1):
            part_id = f"{interval_id}.P{part_number}"
            part_ids.append(part_id)
            part_rows.append(
                {
                    "part_id": part_id,
                    "interval_id": interval_id,
                    "gate1": gate1,
                    "gate2": gate2,
                    "payload_sequence": part.payload_sequence,
                    "ordered_sequence": part.ordered_sequence,
                    "source_type": part.source_type,
                    "is_wt": part.is_wt,
                    "ordered_length": len(part.ordered_sequence),
                    "cost_cents": part.cost_cents,
                    "cost_eur": part.cost_cents / 100,
                    "usage_count": len(part.variant_ids),
                    "primer_length": part.primer_length,
                    "primer_count": part.primer_count,
                }
            )

        for variant_id, part_index in plan.variant_part_indices:
            part = plan.parts[part_index]
            bom_rows.append(
                {
                    "variant_id": variant_id,
                    "interval_id": interval_id,
                    "part_id": part_ids[part_index],
                    "gate1": gate1,
                    "gate2": gate2,
                    "is_wt": part.is_wt,
                    "source_type": part.source_type,
                }
            )

    parts = pd.DataFrame(part_rows)
    bom = pd.DataFrame(bom_rows)
    variant_count = len(variant_dnas)
    non_common_part_ids = set(parts.loc[parts.usage_count < variant_count, "part_id"])
    total_cost_cents = int(parts.cost_cents.sum())
    summary = {
        "requested_variants": variant_count,
        "unique_physical_parts": len(parts),
        "unique_synthesized_non_wt_oligos": int(
            ((parts.source_type == "individual_oligo") & ~parts.is_wt).sum()
        ),
        "directly_synthesized_wt_oligos": int(
            ((parts.source_type == "individual_oligo") & parts.is_wt).sum()
        ),
        "wt_pcr_fragments": int((parts.source_type == "wt_pcr").sum()),
        "terminal_free_parts": int((parts.source_type == "terminal_free").sum()),
        "total_ordered_oligo_nt": int(
            parts.loc[parts.source_type == "individual_oligo", "ordered_length"].sum()
        ),
        "total_primer_nt": int((parts.primer_length * parts.primer_count).sum()),
        "total_cost_cents": total_cost_cents,
        "total_cost_eur": total_cost_cents / 100,
        "total_fragment_usages": len(bom),
        "average_non_wt_parts_per_variant": float((~bom.is_wt).sum()) / variant_count,
        "average_non_common_parts_per_variant": float(
            bom.part_id.isin(non_common_part_ids).sum()
        )
        / variant_count,
        "golden_gate_fragments_per_construct": len(path) - 1,
        "golden_gate_junctions": max(0, len(path) - 2),
        "junctions": [
            {"idx": gate.idx, "bps": gate.bps}
            for gate in path
            if not gate.src_or_target
        ],
    }
    if fidelity_metrics is not None:
        summary["fidelity"] = fidelity_metrics

    return ExactConstructOutput(parts=parts, bom=bom, summary=summary)
