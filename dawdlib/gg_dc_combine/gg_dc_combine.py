import json
from itertools import chain, product
from operator import attrgetter
from typing import Generator, Iterable, Iterator, List, NamedTuple, Tuple, Union

import numpy as np
import pandas as pd

from dawdlib.degenerate_dna.utils import parse_degenerate_codon_csv
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.utils import OligoTableEntry, gate_df_list, parse_dna


class CodonHolder(NamedTuple):
    dna_idx: int
    codons: Iterable[Union[str, None]]

    @property
    def idx(self):
        return self.dna_idx - 1


# def find_oligos(gg_df: pd.DataFrame) -> List[Tuple[int, int]]:
#     """
#
#     Args:
#         gg_df:
#
#     Returns:
#
#     """
#     oligos: List[Tuple[int, int]] = []
#     for ind1, ind2 in zip(gg_df.index[:-1], gg_df.index[1:]):
#         oligos.append((gg_df.loc[ind1, "idx"], gg_df.loc[ind2, "idx"] - 1))
#     return oligos
#
#
# def find_codons_for_oligo(
#     oligo: Tuple[int, int], dc_df: pd.DataFrame
# ) -> List[List[Tuple[int, str]]]:
#     """
#     returns all required combinations of the degenerate codons
#     Args:
#         oligo (Tuple[int, int]): beginning and end positions of the required oligo
#         dc_df (str): degenerate codon data frame
#
#     Returns:
#         a list of lists, each of which has all the required degenerate codons as position and codon tuples.
#     """
#     oligo_codons: List[List[Tuple[int, str]]] = []
#     sub_df = dc_df.loc[
#         dc_df["DNA_POS"].between(oligo[0], oligo[1]),
#         ["DNA_POS"] + [c for c in dc_df.columns if c.startswith("AMBIGUOUS_CODONS")],
#     ]
#     pos_codon: List[List[Tuple[int, str]]] = []
#     for _, row in sub_df.iterrows():
#         pos_codon.append([])
#         for col in sub_df.columns[1:]:
#             if not row[col] or row[col] is np.nan:
#                 continue
#             pos_codon[-1].append((row["DNA_POS"], row[col]))
#
#     for prod in product(*pos_codon):
#         combination: List[Tuple[int, str]] = []
#         for (pos, codon) in prod:
#             combination.append((pos, codon))
#         oligo_codons.append(combination)
#     return oligo_codons
#
#
# def create_dc_oligo(
#     dna: str, pos_codons: List[Tuple[int, str]], oligo: Tuple[int, int]
# ) -> str:
#     """
#     create the DNA string with the required codons
#     Args:
#         dna (str): the full gene DNA
#         pos_codons: list of Gate tuples describing the solution
#         oligo: gate positions for the beginning and end of this oligo
#
#     Returns:
#         (str): an oligo with the degenerate codons
#
#     """
#     dna_copy = dna
#     for (pos, codon) in pos_codons:
#         dna_copy = dna_copy[: pos - 1] + codon + dna_copy[pos + 3 - 1 :]
#     return dna_copy[oligo[0] : oligo[1] + 4]
#
#
# def create_to_order_df(
#     gate_path: List[Gate], deg_df: pd.DataFrame, dna: str, prefix: str, suffix: str
# ) -> pd.DataFrame:
#     """
#     combine the Gate list and degenerate codon table to design the DNA oligos required.
#     Args:
#         gate_path:  list of Gates describing the best path
#         deg_df:  dataframe of the degenerate codons
#         dna:
#         prefix:
#         suffix:
#
#     Returns:
#
#     """
#     oligo_entries = []
#     for gate1, gate2 in zip(gate_path[1:-2], gate_path[2:-1]):
#         oligo_codons = find_codons_for_oligo((gate1.idx, gate2.idx), deg_df)
#         wt_dna = create_dc_oligo(dna, [], (gate1.idx, gate2.idx))
#         wt_entry = OligoTableEntry(
#             wt=True,
#             const=oligo_codons == [[]],
#             full_oligo_dna=prefix + wt_dna + suffix,
#             gate1=gate1,
#             gate2=gate2,
#             gate_gate_dist=gate2.idx - gate1.idx + 3,
#             name=f"{gate1.idx}-{gate2.idx}.wt_dna",
#             oligo_codons=[],
#             oligo_dna=wt_dna,
#         )
#         oligo_entries.append(wt_entry)
#         for ind, oligo_codon in enumerate(oligo_codons):
#             oligo_dna = create_dc_oligo(dna, oligo_codon, (gate1.idx, gate2.idx))
#             if not oligo_codon:
#                 continue
#             oligo_entries.append(
#                 wt_entry._replace(
#                     wt=False,
#                     oligo_dna=oligo_dna,
#                     full_oligo_dna=prefix + oligo_dna + suffix,
#                     name=f"{gate1.idx}-{gate2.idx}.{len(oligo_codons)}.{ind}",
#                     oligo_codons=oligo_codon,
#                 )
#             )
#     return pd.DataFrame.from_records(oligo_entries, columns=OligoTableEntry._fields)


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
    g2_idx = g2.span()[1]
    try:
        g2_idx += 1
    except TypeError:
        pass
    o_slice = slice(g1.idx, g2_idx)
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
    idxs = [cdn.idx for cdn in cdns]
    muts_cdns = [cdn.codons for cdn in cdns]
    for i, m_cdns in enumerate(product(*muts_cdns)):
        e_dna = ""
        e_cdns = []
        for idx1, idx2, cdn in zip(idxs[:-1], idxs[1:], m_cdns[:-1]):
            if cdn is None:
                e_dna += dna[idx1:idx2]
            else:
                e_cdns.append(cdn)
                e_dna += cdn + dna[idx1 + 3 : idx2]
        yield wt._replace(
            wt=False,
            oligo_dna=e_dna,
            name=f"{idxs[0]}-{idxs[-1]}.{i+1}",
            oligo_codons=e_cdns,
            const=False,
        )


def find_cdns(
    g1: Gate, g2: Gate, dna: str, deg_codons: Iterable[CodonHolder],
) -> List[CodonHolder]:
    """
    Given indices returns the codons which are between the indices
    Args:
        indices (Tuple[int, int, int]): A tuple of [start, stop, step] (produced by slice().indices())
        deg_codons (Iterable[CodonHolder]): the entire codons

    Returns:
        List[CodonHolder]: Matching codons between start and stop

    """
    o_slice = slice(g1.idx, g2.span()[1])
    indices = o_slice.indices(len(dna))
    return [x for x in deg_codons if indices[0] < x.idx < indices[1]]


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
    g1_idx = g1.idx
    try:
        g1_idx += 1
    except TypeError:
        pass
    g2_idx = g2.span()[1]
    try:
        g2_idx += 2
    except TypeError:
        pass
    g1_edge_codon = CodonHolder(g1_idx, [None])
    g2_edge_codon = CodonHolder(g2_idx, [None])
    cdns += [g1_edge_codon, g2_edge_codon]
    return sorted(cdns, key=attrgetter("idx"))


def create_oligo(
    g1: Gate, g2: Gate, cdns: List[CodonHolder], dna: str, prefix: str, suffix: str
) -> Iterable[OligoTableEntry]:
    oligos = []
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
        o_pre = prefix if g1.idx is not None else ""
        o_suf = suffix if g2.idx is not None else ""
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
