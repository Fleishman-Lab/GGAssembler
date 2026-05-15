import json
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple

import pandas as pd
from Bio.Seq import Seq

from dawdlib.degenerate_dna.codon_selector import CodonSelector
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.utils import check_for_restriction_sites, find_dna_var_poss
from dawdlib.gg_dc_combine.gg_dc_combine import get_wt


class Mutation(NamedTuple):
    aa_pos: int
    from_aa: str
    to_aa: str


class Variant(NamedTuple):
    variant_id: str
    target_aas: str
    mutations: Sequence[Mutation]
    total_score: Optional[float] = None


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
    result = {}

    for variant in variants:
        dna = list(wt_dna)
        for mutation in variant.mutations:
            dna_idx = (mutation.aa_pos - 1) * 3
            codon = choose_exact_codon(codon_selector, mutation.to_aa)
            dna[dna_idx : dna_idx + 3] = codon

        variant_dna = "".join(dna)
        if enzymes:
            ok, enzyme, sites = check_for_restriction_sites(variant_dna, list(enzymes))
            if not ok:
                raise ValueError(
                    f"{variant.variant_id} introduces {enzyme} sites at {sites}"
                )
        result[variant.variant_id] = variant_dna

    return result


def exact_variant_oligo_table(
    path: Sequence[Gate],
    wt_dna: str,
    variant_dnas: Dict[str, str],
    prefix: str,
    suffix: str,
    project_name: str = "",
) -> pd.DataFrame:
    rows = []
    cassette_id = 0
    name_prefix = f"{project_name}." if project_name else ""

    for g1, g2 in zip(path[:-1], path[1:]):
        wt_entry = get_wt(g1, g2, const=True, dna=wt_dna)
        o_pre = prefix if not g1.src_or_target else ""
        o_suf = suffix if not g2.src_or_target else ""
        wt_segment = wt_entry.oligo_dna
        variant_segments = {
            variant_id: get_wt(g1, g2, const=True, dna=variant_dna).oligo_dna
            for variant_id, variant_dna in variant_dnas.items()
        }
        is_variable = any(segment != wt_segment for segment in variant_segments.values())
        current_cassette_id = None
        if is_variable:
            cassette_id += 1
            current_cassette_id = cassette_id

        for variant_id, segment in variant_segments.items():
            entry = wt_entry._replace(
                wt=segment == wt_segment,
                const=not is_variable,
                oligo_dna=segment,
                full_oligo_dna=o_pre + segment + o_suf,
                name=f"{name_prefix}{variant_id}.{g1.idx}-{g2.idx}",
            )
            rows.append(
                {
                    **entry._asdict(),
                    "variant_id": variant_id,
                    "cassette_id": current_cassette_id,
                }
            )

    df = pd.DataFrame(rows)
    df["gate1"] = df["gate1"].map(json.dumps)
    df["gate2"] = df["gate2"].map(json.dumps)
    return df
