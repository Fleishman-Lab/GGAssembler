import json
from itertools import chain, product
from operator import attrgetter
from typing import Generator, Iterable, Iterator, List, NamedTuple, Optional, Union

import pandas as pd

from dawdlib.degenerate_dna.utils import parse_degenerate_codon_csv
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.utils import OligoTableEntry, gate_df_list, parse_dna


class CodonHolder(NamedTuple):
    dna_idx: Optional[int]
    codons: Iterable[Union[str, None]]

    @property
    def idx(self):
        if self.dna_idx is not None:
            return self.dna_idx - 1
        return self.dna_idx


def dc_df_codon_list(dc_df: pd.DataFrame) -> Iterable[CodonHolder]:
    """
    Used to transform a dataframe of degenerate codons with some line having empty
    AMBIGUOUS_CODONS columns into CodonHolders
    Args:
        dc_df (pd.Dataframe): The degenerate codon dataframe

    Returns:
        iterable[CodonHolder]:

    """
    codon_idx = dc_df.columns.to_list().index("AMBIGUOUS_CODONS1")
    return [
        CodonHolder(tup.DNA_POS, list(filter(None, tup[codon_idx:])))
        for tup in dc_df.itertuples(index=False)
    ]


def get_wt(g1: Gate, g2: Gate, const: bool, dna: str,) -> OligoTableEntry:
    """
    Returns the wt oligo between the two gates
    Args:
        g1 (Gate): First gate
        g2 (Gate): Second gate
        const (bool): Is this oligo constant
        dna (str): Full gene dna from which the appropriate section will be taken

    Returns:
        OligoTableEntry: Of the WT sequence

    """
    g1_idx = g1.idx if not g1.src_or_target else None
    g2_idx = g2.span()[1] + 1 if not g2.src_or_target else None
    o_slice = slice(g1_idx, g2_idx)
    indices = o_slice.indices(len(dna))
    wt_entry = OligoTableEntry(
        wt=True,
        const=const,
        full_oligo_dna=dna[o_slice],
        gate1=g1,
        gate2=g2,
        gate_gate_dist=indices[1] - indices[0],
        name=f"{g1.idx}-{g2.idx}.wt_dna",
        oligo_codons=[],
        oligo_dna=dna[o_slice],
    )
    return wt_entry


def yield_mut(
    cdns: List[CodonHolder], dna: str, wt: OligoTableEntry
) -> Generator[OligoTableEntry, None, None]:
    """
    Yields mutant oligos which is a combination of codons provided in cdns.
    When a codon is None the wt sequence is used otherwise the codon string
    replaces the right DNA segment.
    Args:
        cdns (List[CodonHolder]): List of codons for the right oligo which
                                  includes at the start and end the oligo
                                  gates with a None codon.
        dna (str): DNA of the entire gene
        wt (OligoTableEntry): The WT oligo of the current segment (from gate1 to gate 2)

    Yields:
        OligoTableEntry: The appropriate oligo
    """
    muts_cdns = [cdn.codons for cdn in cdns]
    for i, m_cdns in enumerate(product(*muts_cdns)):
        e_dna = ""
        e_cdns = []
        for ch1, ch2, cdn in zip(cdns[:-1], cdns[1:], m_cdns[:-1]):
            dna_slice = slice(ch1.idx, ch2.idx)
            idx1, idx2, _ = dna_slice.indices(len(dna))
            if cdn is None:
                e_dna += dna[idx1:idx2]
            else:
                e_cdns.append(cdn)
                e_dna += cdn + dna[idx1 + 3 : idx2]
        yield wt._replace(
            wt=False,
            oligo_dna=e_dna,
            name=f"{cdns[0].idx}-{cdns[-1].idx}.{i+1}",
            oligo_codons=e_cdns,
            const=False,
        )


def find_cdns(
    g1: Gate, g2: Gate, dna: str, deg_codons: Iterable[CodonHolder],
) -> List[CodonHolder]:
    """
    Given indices returns the codons which are between the indices
    Args:
        g1 (Gate): First gate
        g2 (Gate): Second gate
        dna (str): The gene DNA string
        deg_codons (Iterable[CodonHolder]): the entire codons

    Returns:
        List[CodonHolder]: Matching codons between start and stop

    """
    g1_idx = g1.idx if not g1.src_or_target else None
    g2_idx = g2.span()[1] if not g2.src_or_target else None
    o_slice = slice(g1_idx, g2_idx)
    start, end, _ = o_slice.indices(len(dna))
    if g1.src_or_target:
        start -= 1
    if g2.src_or_target:
        end += 1
    return [x for x in deg_codons if start < x.idx < end]


def cdn_add_gates(g1: Gate, g2: Gate, cdns: List[CodonHolder]) -> List[CodonHolder]:
    """
    Add gates g1 and g2 to the codon list and sort them by index(idx)
    Args:
        g1 (Gate): First gate of segment (start)
        g2 (Gate): Second gate of segment (stop)
        cdns (List[CodonHolder])): codon list to add gates to

    Returns:
        List[CodonHolder]: A sorted list of codons by index with gate1 and gate2 added to it

    """
    cdns = sorted(cdns, key=attrgetter("idx"))
    start = CodonHolder(g1.idx + 1, [None])
    if g1.src_or_target:
        start = start._replace(dna_idx=None)
    cdns.insert(0, start)
    end = CodonHolder(g2.span()[1] + 2, [None])
    if g2.src_or_target:
        end = end._replace(dna_idx=None)
    cdns.append(end)
    return cdns


def create_oligo(
    g1: Gate, g2: Gate, cdns: List[CodonHolder], dna: str, prefix: str, suffix: str
) -> List[Union[List[OligoTableEntry], Iterable[OligoTableEntry]]]:
    oligos: List[Union[List[OligoTableEntry], Iterable[OligoTableEntry]]] = []
    const = 0 == len(cdns)
    wt = get_wt(g1, g2, const, dna)
    oligos.append([wt._replace(full_oligo_dna=prefix + wt.oligo_dna + suffix,)])
    if not const:
        cdns = cdn_add_gates(g1, g2, cdns)
        oligos.append(
            map(
                lambda x: x._replace(full_oligo_dna=prefix + x.oligo_dna + suffix,),
                yield_mut(cdns, dna, wt),
            )
        )
    return oligos


def create_oligos(
    pth: List[Gate],
    deg_codons: Iterable[CodonHolder],
    dna: str,
    prefix: str,
    suffix: str,
    project_name: str = "",
) -> Iterator[OligoTableEntry]:
    """
    Creates all oligo entries for a given path, degenerate codons and dna
    Args:
        pth (List[Gate]): Specifies the gates, the location of oligos start and end.
        deg_codons (Iterable[CodonHolder]): Specifies all of the degenerate codons
                                            needed to be embedded into the DNA.
        dna (str): The entire gene DNA
        prefix (str): The DNA prefix to add to every oligo produced
        suffix (suffix): The DNA suffix to add to every oligo produced
        project_name (str): A name to be added as a prefix to oligo name

    Returns:
        Iterator[OligoTableEntry]: An Iterator over oligo entries.

    """
    oligos = []
    if project_name:
        project_name += "."
    for g1, g2 in zip(pth[:-1], pth[1:]):
        cdns = find_cdns(g1, g2, dna, deg_codons)
        o_pre = prefix if not g1.src_or_target else ""
        o_suf = suffix if not g2.src_or_target else ""
        oligos.extend(create_oligo(g1, g2, cdns, dna, o_pre, o_suf))
    return map(
        lambda x: x._replace(name=f"{project_name}{x.name}"),
        chain.from_iterable(oligos),
    )


def gate_cdn_oligos(
    pth: List[Gate],
    deg_codons: Iterable[CodonHolder],
    dna: str,
    prefix: str,
    suffix: str,
    project_name: str = "",
) -> pd.DataFrame:
    """
    A helper function create a DataFrame based on oligo entries and to parse correctly gates into json.
    See Also:
        create_oligos
    Returns:
        pd.DataFrame: whose entries are OligoEntry from create_oligos

    """
    to_order_df = pd.DataFrame.from_records(
        create_oligos(pth, deg_codons, dna, prefix, suffix, project_name),
        columns=OligoTableEntry._fields,
    )
    to_order_df["gate1"] = to_order_df["gate1"].map(json.dumps)
    to_order_df["gate2"] = to_order_df["gate2"].map(json.dumps)
    return to_order_df


def combine_gate_path_deg_codons(
    dc_table_file: str,
    gate_path_file: str,
    dna_file: str,
    prefix: str,
    suffix: str,
    to_order_df_file: str,
    project_name: str = "",
) -> None:
    """
    A helper function which loads the right files and uses gate_cdn_oligos to write the dataframe into a file.
    Args:
        dc_table_file (str): A path pointing to a degenerate table csv
        gate_path_file (str: A path pointing to a gate path table csv
        dna_file (str): A path pointing to dna in fasta format
        prefix (str): The prefix to add to every oligo
        suffix (str): The prefix to add to every oligo
        to_order_df_file (str):  A path pointing to the desired oligo tale output
        project_name (str): A name to be added as a prefix to oligo name
    """
    dna = parse_dna(dna_file)
    deg_codons = dc_df_codon_list(parse_degenerate_codon_csv(dc_table_file))
    gate_path = gate_df_list(pd.read_csv(gate_path_file))
    to_order_df = gate_cdn_oligos(
        gate_path, deg_codons, dna, prefix, suffix, project_name
    )
    to_order_df.to_csv(to_order_df_file)
