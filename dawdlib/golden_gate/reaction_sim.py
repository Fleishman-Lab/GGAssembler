from itertools import chain, combinations, product
from typing import Generator, Iterable, Iterator, List, NamedTuple, Tuple, Union

import networkx as nx
import numpy as np
import pandas as pd
from Bio.Restriction import Restriction
from Bio.Restriction.Restriction import Ov3, Ov5
from Bio.Seq import Seq
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.utils import OligoTableEntry, Requirements


class OverHang(NamedTuple):
    """
    A NamedTuple representing the length of hanging dna sections.
    """

    fprime: int = 0
    tprime: int = 0


class OverHangDNA(NamedTuple):
    """
    A NamedTuple representing the DNA strings of hanging dna sections.
    """

    fprime: str = ""
    tprime: str = ""


class SDNASection(NamedTuple):
    """
    A NamedTuple representing a single strand of DNA
    """

    dna: str = ""
    start: int = 0
    end: int = 0
    hang: OverHang = OverHang()
    hang_dna: OverHangDNA = OverHangDNA()


class DDNASection(NamedTuple):
    """
    A NamedTuple representing a double strand of DNA
    """

    fprime: SDNASection = SDNASection()
    tprime: SDNASection = SDNASection()
    is_wt: bool = False

    def get_hang_dna(self) -> Iterator[str]:
        """

        Returns:
            Iterator[str]: Iterator over the hanging DNA.
        """
        return chain(self.fprime.hang_dna, self.tprime.hang_dna)


class ReactionGraph(nx.Graph):
    """
    A wrapper class around ~networkx.Graph which has a source and target fields.
    """

    source: DDNASection = DDNASection(is_wt=True)
    target: DDNASection = DDNASection(
        fprime=SDNASection(start=np.iinfo(np.int).max, end=np.iinfo(np.int).max),
        is_wt=True,
    )


class ReactionSim:

    reaction_graph: ReactionGraph

    def __init__(self, gg_data: GGData, reqs: Requirements, enzymes: List[str]):
        self.gg_data = gg_data
        self.reqs = reqs
        self.enzymes: List[Union[Ov3, Ov5]] = []
        self.wt_oligos: List[DDNASection]
        for enzyme in enzymes:
            enz_cls = getattr(Restriction, enzyme, None)
            if enz_cls is not None:
                self.enzymes.append(enz_cls)

    def load_oligo_table(self, table_path: str) -> Iterator[DDNASection]:
        """
        Turns a table of ~dawdlib.goldengate.utils.OligoTableEntry into corresponding ~DDNASection objects
        Args:
            table_path (str): The path to a CSV file of OligoTableEntries.

        Returns:
            Iterator[DDNASection]: An iterator over corresponding ~DDNASection objects

        """
        oligo_table = pd.read_csv(table_path)
        ddna_nodes: List[Generator[DDNASection, None, None]] = []
        for enzyme in self.enzymes:
            for row in oligo_table.itertuples(index=False):
                oligo_entry = OligoTableEntry(**row._asdict())
                dna_segs_it = enzyme_dna_segments(enzyme, oligo_entry.full_oligo_dna)
                ddna_nodes.append(
                    get_ddna_sections(
                        enzyme, dna_segs_it, oligo_entry.gate1.idx, oligo_entry.wt
                    )
                )
        return chain.from_iterable(ddna_nodes)

    def create_reaction_graph(self, table_path: str) -> None:
        """
        Creates the reaction graph given a table of ~dawdlib.goldengate.utils.OligoTableEntry
        Args:
            table_path (str): The path to a CSV file of OligoTableEntries.
        """
        ddna_iter = self.load_oligo_table(table_path)
        rgraph = ReactionGraph()
        rgraph.add_nodes_from(ddna_iter)
        rgraph.add_edges_from(
            get_compatible_oligos(
                self.gg_data, iter(rgraph.nodes), self.reqs.gate_crosstalk_max
            )
        )
        sources, targets = find_sources_targets(iter(rgraph.nodes))
        rgraph.add_edges_from(
            gen_src_trgt_edges(rgraph.source, rgraph.target, sources, targets)
        )
        self.reaction_graph = rgraph

    def get_wt_dna(self) -> Generator[str, None, None]:
        """
        Finds all WT dna segments that can be connected in the reaction.

        Yields:
            str: The WT DNA sequence found in the reaction graph

        """
        sub_reaction_g = self.reaction_graph.subgraph(
            filter(lambda x: x.is_wt, self.reaction_graph.nodes)
        ).copy()
        pth: List[DDNASection]
        for pth in nx.all_simple_paths(
            sub_reaction_g, self.reaction_graph.source, self.reaction_graph.target
        ):
            yield "".join((node.fprime.dna for node in pth))

    def verify_reaction(self, product_len: int = 0) -> bool:
        """
        Checks if all DNA sequences in the reaction are connected correctly
        by verifying the numbering of the DNA product is in ascending order.

        If product_len is provided also verifies that the DNA produced
        has length exactly equal to product_len.

        Args:
            product_len (int): Allows to verify all products have len of exactly product_len base pairs.

        Returns:
            bool: True

        """
        pth: List[DDNASection]
        for pth in nx.all_simple_paths(
            self.reaction_graph, self.reaction_graph.source, self.reaction_graph.target
        ):
            if not all(
                (
                    pth[i].fprime.start < pth[i + 1].fprime.start
                    and pth[i].fprime.end < pth[i].fprime.end
                    for i in range(len(pth) - 1)
                )
            ):
                return False
            if 0 < product_len != len("".join(node.fprime.dna for node in pth)):
                return False
        return True


def gen_src_trgt_edges(
    source: DDNASection,
    target: DDNASection,
    sources: Iterable[DDNASection],
    targets: Iterable[DDNASection],
) -> Iterator[Tuple[DDNASection, DDNASection]]:
    """
    Generates the required edges so the graph structure will have a single source and target
    Args:
        source (DDNASection): source node
        target (DDNASection): target node
        sources (List[DDNASection]): current sources
        targets (List[DDNASection]): current targets

    Returns:
        Iterator[Tuple[DDNASection, DDNASection]]: An iterator over tuples of pairs for the new source and target

    """
    return chain(product([[source], sources]), product([[target]], targets))


def find_sources_targets(
    ddnas: Iterator[DDNASection],
) -> Tuple[Iterable[DDNASection], Iterable[DDNASection]]:
    """
    Finds all sources and targets of the given DDNASections.
    Sources are sections with minimal fprime.start targets are sections with maximal fprime.end

    Args:
        ddnas (Iterator[DDNASection]): An iterator over the DDNASection to search over.

    Returns:
        Tuple[Iterable[DDNASection], Iterable[DDNASection]]: A tuple of the sources, targets.

    """
    min_source_start = np.Inf
    max_target_end = np.NINF
    sources: List[DDNASection] = []
    targets: List[DDNASection] = []
    for ddna in ddnas:
        if ddna.fprime.start < min_source_start:
            min_source_start = ddna.fprime.start
            sources = [ddna]
        elif ddna.fprime.start == min_source_start:
            sources.append(ddna)
        if ddna.fprime.end > max_target_end:
            max_target_end = ddna.fprime.end
            targets = [ddna]
        elif ddna.fprime.end == max_target_end:
            targets.append(ddna)
    return sources, targets


def get_compatible_oligos(
    gg_data: GGData, ddnas: Iterator[DDNASection], gate_crosstalk_max: int = 1000
) -> Generator[Tuple[DDNASection, DDNASection], None, None]:
    """

    Args:
        gg_data (GGData): A wrapper around the gates data
        ddnas (Iterator[DDNASection]):
        gate_crosstalk_max (int): gg_data.gates_all_scores above this score marks the two dna sections as compatible

    Yields:
         Tuple[DDNASection, DDNASection] - Signifying these two dna sections can connect to one another in the reaction
    """
    dtup: Tuple[DDNASection, ...]
    for dtup in combinations(ddnas, 2):
        d1, d2 = dtup
        for h1, h2 in product(d1.get_hang_dna(), d2.get_hang_dna()):
            if gg_data.gates_all_scores(h1, h2) >= gate_crosstalk_max:
                yield d1, d2
                continue


def get_ddna_sections(
    enzyme: Union[Ov5, Ov3],
    dna_segs: Iterator[Tuple[str, str]],
    start: int,
    is_wt: bool = False,
) -> Generator[DDNASection, None, None]:
    """

    Args:
        enzyme (Ov5): Restriction enzyme which was used to digest the DNA.
        dna_segs (Iterator[Tuple[str, str]]): An iterator over double stranded DNA segments
                                              (after digestion by restrction enzyme).
        start (int): An index pointing to the start of the DNA segment or the location of the gate.
        is_wt (bool): Whether the current oligo is WT or not.

    Yields:
        DDNASection: The double stranded DNA segment.
    """
    prev_ddna = DDNASection(fprime=SDNASection(end=start))
    for fe, te in dna_segs:
        f_overhang, t_overhang = find_overhang(
            len(fe), len(te), prev_ddna.fprime.hang.tprime, prev_ddna.tprime.hang.fprime
        )
        f_dna_overhang, t_dna_overhang = find_overhand_dna(
            fe, te, f_overhang, t_overhang
        )
        prev_ddna = DDNASection(
            fprime=SDNASection(
                fe,
                prev_ddna.fprime.end,
                prev_ddna.fprime.end + len(fe),
                f_overhang,
                f_dna_overhang,
            ),
            tprime=SDNASection(dna=te, hang=t_overhang, hang_dna=t_dna_overhang),
            is_wt=is_wt,
        )
        # Skip yielding sections which have the enzyme recognition site
        if enzyme.search(fe) or enzyme.search(te):
            continue
        yield prev_ddna


def enzyme_dna_segments(enzyme: Union[Ov5, Ov3], dna: str) -> Iterator[Tuple[str, str]]:
    """

    Args:
        enzyme (Ov5): Restriction enzyme used to cut the dna
        dna (str): 5' DNA string

    Returns:
        Iterator[Tuple[str, str]]: An iterator over the double stranded DNA segments
    """
    seq = Seq(dna)
    f2t = map(str, enzyme.catalyse(seq))
    t2f = reversed(
        list(
            map(
                lambda x: "".join(reversed(x)),
                enzyme.catalyse(seq.reverse_complement()),
            )
        )
    )
    return iter(zip(f2t, t2f))


def find_overhang(
    flen: int, tlen: int, pft: int, ptf: int
) -> Tuple[OverHang, OverHang]:
    """
    Find DNA hanging ends.

    Takes care of six different cases:
        1. ----------
           ---------------
           flen: 0, tlen: 0 -> flen: 0, tlen: 5
        2. ---------------
                ----------
           flen: 0, tlen: 5 -> flen: 0, tlen: 0
        3. ---------------
                ---------------
           flen: 0, tlen: 5 -> flen: 0, tlen: 5
        And the reverse of these cases
    Args:
        flen (int): Length of DNA string of 5' to 3'
        tlen (int): Length of DNA string of 3' to 5'
        pft (int): previous 5' segment 3' side hanging length
        ptf (int): previous 3' segment 5' side hanging length

    Returns:
        Tuple[OverHang, OverHang]: A tuple with the lengths of the hanging DNA of both stands (5' and 3').
    """

    cft = flen + pft - tlen
    cft = cft if cft >= 0 else 0
    ctf = tlen + ptf - flen
    ctf = ctf if ctf >= 0 else 0
    return OverHang(ptf, cft), OverHang(pft, ctf)


def find_overhand_dna(
    fdna: str, tdna: str, f_overhang: OverHang, t_overhang: OverHang
) -> Tuple[OverHangDNA, OverHangDNA]:
    """
    See ~find_overhang for cases and logic.

    Args:
        fdna (str): The 5' DNA string
        tdna (str): The 3' DNA string:
        f_overhang (OverHang): The 5' strand overhangs
        t_overhang (OverHang): The 5' strand overhangs

    Returns:
        Tuple[OverHangDNA, OverHangDNA]
    """
    return (
        OverHangDNA(fdna[: t_overhang.tprime], fdna[: -f_overhang.tprime]),
        OverHangDNA(tdna[: f_overhang.fprime], tdna[: -t_overhang.fprime]),
    )
