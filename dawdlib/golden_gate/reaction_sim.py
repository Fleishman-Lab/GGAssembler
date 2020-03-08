import json
from itertools import chain, combinations, product
from typing import Dict, Generator, Iterable, Iterator, List, NamedTuple, Tuple, Union

import networkx as nx
import numpy as np
import pandas as pd
from Bio.Restriction import Restriction
from Bio.Restriction.Restriction import Ov3, Ov5
from Bio.Seq import Seq
from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.utils import OligoTableEntry, Requirements, ambiguous_dna_unambiguous
from tabulate import tabulate


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

    # 5' to 3'
    fwd: SDNASection = SDNASection()
    rev: SDNASection = SDNASection()
    is_wt: bool = False

    def get_hang_dna(self) -> Iterator[str]:
        """

        Returns:
            Iterator[str]: Iterator over the hanging DNA.
        """
        return iter(
            (
                self.fwd.hang_dna.fprime,
                self.fwd.hang_dna.tprime,
                self.rev.hang_dna.tprime,
                self.rev.hang_dna.fprime,
            )
        )

    def get_hang_dna_rev(self) -> Iterator[str]:
        """

        Returns:
            Iterator[str]: Iterator over the hanging DNA.
            The reverse strand is reversed to match supplied tables.
        """
        return iter(
            (
                self.fwd.hang_dna.fprime,
                self.fwd.hang_dna.tprime,
                self.rev.hang_dna.tprime[::-1],
                self.rev.hang_dna.fprime[::-1],
            )
        )


class ReactionGraph(nx.Graph):
    """
    A wrapper class around ~networkx.Graph which has a source and target fields.
    """

    source: DDNASection = DDNASection(is_wt=True)
    target: DDNASection = DDNASection(
        fwd=SDNASection(start=np.iinfo(np.int).max, end=np.iinfo(np.int).max),
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
        oligo_table = pd.read_csv(
            table_path,
            header=0,
            names=OligoTableEntry._fields,
            converters={
                "gate1": lambda x: Gate(*json.loads(x)),
                "gate2": lambda x: Gate(*json.loads(x)),
            },
        )
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

    def get_wt_dna(self) -> Generator[Tuple[str, int, int, int], None, None]:
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
            dna = "".join((node.fwd.dna for node in pth))
            start = 0
            end = 0
            try:
                start = pth[pth.index(self.reaction_graph.source) + 1].fwd.start
                end = pth[pth.index(self.reaction_graph.target) - 1].fwd.end
            except ValueError:
                pass
            except IndexError:
                pass
            yield dna, len(pth), start, end

    def verify_reaction(
        self, expected_prod_len: int = 0, expected_no_segments: int = 0
    ) -> Tuple[bool, Union[int, List[DDNASection]]]:
        """
        Checks if all DNA sequences in the reaction are connected correctly
        by verifying the numbering of the DNA product is in ascending order.

        If product_len is provided also verifies that the DNA produced
        has length exactly equal to product_len.

        Args:
            expected_prod_len (int): Allows to verify all products have len of exactly expected_prod_len base pairs.
            expected_no_segments (int): Allows to verify the product is composed of exactly
                expected_no_segments segments.

        Returns:
            Tuple: 1. in case all products are correctly assembled (True, number of products)
                   2. Otherwise: (False, the product which is not correct

        """
        pth: List[DDNASection]
        cnt = -1
        for cnt, pth in enumerate(
            nx.all_simple_paths(
                self.reaction_graph,
                self.reaction_graph.source,
                self.reaction_graph.target,
            )
        ):
            if not all(
                (
                    pth[i].fwd.start < pth[i + 1].fwd.start
                    and pth[i].fwd.end < pth[i + 1].fwd.end
                    for i in range(len(pth) - 1)
                )
            ):
                return False, pth
            if 0 < expected_prod_len != len("".join(node.fwd.dna for node in pth)):
                return False, pth
            if 0 < expected_no_segments != len(pth):
                return False, pth
        return True, cnt + 1

    def get_all_products(self) -> Generator[Tuple[str, int, int, int], None, None]:
        """
        Use wisely! this method yields all possible products of the golden gate reaction.
        Yields:
            Tuple: [dna, number of segments, start position, end position] 
        """
        pth: List[DDNASection]
        for pth in nx.all_simple_paths(
            self.reaction_graph, self.reaction_graph.source, self.reaction_graph.target
        ):
            yield "".join(node.fwd.dna for node in pth), len(pth), pth[
                0
            ].fwd.start, pth[-1].fwd.end


class ReactionCLI:
    """
    Golden gate reaction simulator.
    Provides two methods to verify golden gate products:
        1. get_wt: Which tries to retrieve the WT sequence.
        2. verify: Which tries to verify the entire reaction.
    """

    sim: ReactionSim

    def _setup(
        self,
        table_path: str,
        enzymes: List[str],
        gate_crosstalk_max: int = 1000,
        neb_table_temp: int = 25,
        neb_table_time: int = 18,
    ):

        self.sim = ReactionSim(
            GGData(neb_table_temp=neb_table_temp, neb_table_time=neb_table_time),
            Requirements(0, 0, 0, gate_crosstalk_max=gate_crosstalk_max),
            enzymes,
        )
        self.sim.create_reaction_graph(table_path)
        return self

    def get_wt(
        self,
        table_path: str,
        enzymes: Union[List[str], str],
        gate_crosstalk_max: int = 1000,
        neb_table_temp: int = 25,
        neb_table_time: int = 18,
    ):
        """
        Prints the WT sequence(s) created by "simulating" the golden gate reaction.

        Args:
            table_path (str): Path to segments table
            enzymes (str|list[str]): which enzyme(s) to use in the reaction given. Enzymes should be given in a list format, i.e. [BsaI]. You must respect the usual naming convention with the upper case letters and Latin numbering (in upper case as well).
            gate_crosstalk_max (int): above which threshold two gates are considered as able to connect (default=1000)
            neb_table_temp (int): Specifies which NEB data to use (default=25)
            neb_table_time (int): Specifies which NEB data to use (default=18)
        """
        if isinstance(enzymes, str):
            enzymes = [enzymes]
        self._setup(
            table_path, enzymes, gate_crosstalk_max, neb_table_temp, neb_table_time
        )
        print("Found WT DNA sequences:")
        for dna, no_segs, dna_start_pos, dna_end_pos in self.sim.get_wt_dna():
            message = (
                f"DNA sequence: \t {dna}\n"
                f"Is composed of {no_segs} segments and spans dna positions {dna_start_pos} to {dna_end_pos}"
            )
            print(message)

    def verify(
        self,
        table_path: str,
        enzymes: Union[List[str], str],
        expected_dna_len: int = 0,
        expected_no_segments: int = 0,
        gate_crosstalk_max: int = 1000,
        neb_table_temp: int = 25,
        neb_table_time: int = 18,
    ):
        """
        Verifies the "correctness" of the products in golden gate simulation by:
            1. Check that oligo number is in ascending order.
                For example: that oligo of BPs 12-60 is before oligo of BPs 61-75.
            2. If expected_dna_len is given,
                verifies all reaction products are of the required length.
            3. If  expected_no_segments is given,
                verifies all reaction products are composed of required number of oligos.
        Args:
            table_path (str): Path to segments table
            enzymes (str|list[str]): which enzyme(s) to use in the reaction given. Enzymes should be given in a list format, i.e. [BsaI]. You must respect the usual naming convention with the upper case letters and Latin numbering (in upper case as well).
            expected_dna_len (int): The number of BPs each product is expected to contain
            expected_no_segments (int): The number of segments each product is expected to contain
            gate_crosstalk_max (int): above which threshold two gates are considered as able to connect (default=1000)
            neb_table_temp (int): Specifies which NEB data to use (default=25)
            neb_table_time (int): Specifies which NEB data to use (default=18)

        """
        if isinstance(enzymes, str):
            enzymes = [enzymes]
        self._setup(
            table_path, enzymes, gate_crosstalk_max, neb_table_temp, neb_table_time
        )
        success, extra = self.sim.verify_reaction(
            expected_prod_len=expected_dna_len,
            expected_no_segments=expected_no_segments,
        )
        if success:
            print(
                f"Reaction simulation is successful, {extra} unique sequences were found."
            )
        else:
            print(
                "The reaction simulation failed the following illegal DNA sequence was found:\n"
            )
            table = [
                list(
                    chain(
                        [dd_seg.fwd.dna, dd_seg.fwd.start, dd_seg.fwd.end],
                        dd_seg.get_hang_dna(),
                    )
                )
                for dd_seg in extra
            ]
            print(
                tabulate(
                    table,
                    headers=[
                        "DNA",
                        "start",
                        "end",
                        "fwd 5' overhang dna",
                        "fwd 3' overhang dna",
                        "rev 3' overhang dna",
                        "rev 5' overhang dna",
                    ],
                )
            )

    def get_all_products(
        self,
        table_path: str,
        enzymes: Union[List[str], str],
        gate_crosstalk_max: int = 1000,
        neb_table_temp: int = 25,
        neb_table_time: int = 18,
    ):
        """
        Use with caution!!!
        This method prints all golden gate reaction products in a table format.
        Args:
            table_path (str): Path to segments table
            enzymes (str|list[str]): which enzyme(s) to use in the reaction given. Enzymes should be given in a list format, i.e. [BsaI]. You must respect the usual naming convention with the upper case letters and Latin numbering (in upper case as well).
            gate_crosstalk_max (int): above which threshold two gates are considered as able to connect (default=1000)
            neb_table_temp (int): Specifies which NEB data to use (default=25)
            neb_table_time (int): Specifies which NEB data to use (default=18)
        """
        if isinstance(enzymes, str):
            enzymes = [enzymes]
        self._setup(
            table_path, enzymes, gate_crosstalk_max, neb_table_temp, neb_table_time
        )
        products = self.sim.get_all_products()
        print(tabulate(products, headers=["DNA", "No. segments", "Start", "End"]))


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
    return chain(product([source], sources), product([target], targets))


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
        if ddna.fwd.start < min_source_start:
            min_source_start = ddna.fwd.start
            sources = [ddna]
        elif ddna.fwd.start == min_source_start:
            sources.append(ddna)
        if ddna.fwd.end > max_target_end:
            max_target_end = ddna.fwd.end
            targets = [ddna]
        elif ddna.fwd.end == max_target_end:
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
        for h1, h2 in product(
            filter(len, d1.get_hang_dna_rev()), filter(len, d2.get_hang_dna_rev()),
        ):
            if gg_data.gates_scores(h1, h2) >= gate_crosstalk_max:
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
    prev_ddna = DDNASection(fwd=SDNASection(end=start))
    for fwd, rev in dna_segs:
        f_overhang, t_overhang = find_overhang(
            len(fwd), len(rev), prev_ddna.fwd.hang, prev_ddna.rev.hang
        )
        f_dna_overhang, t_dna_overhang = find_overhand_dna(
            fwd, rev, f_overhang, t_overhang
        )
        prev_ddna = DDNASection(
            fwd=SDNASection(
                fwd,
                prev_ddna.fwd.end,
                prev_ddna.fwd.end + len(fwd),
                f_overhang,
                f_dna_overhang,
            ),
            rev=SDNASection(dna=rev, hang=t_overhang, hang_dna=t_dna_overhang),
            is_wt=is_wt,
        )
        # Skip yielding sections which have the enzyme recognition site
        if enzyme.compsite.search(fwd) is not None:
            prev_ddna = prev_ddna._replace(fwd=prev_ddna.fwd._replace(end=start))
            continue
        yield prev_ddna


def enzyme_dna_segments(
    enzyme: Union[Ov5, Ov3], udna: str
) -> Iterator[Tuple[str, str]]:
    """

    Args:
        enzyme (Ov5): Restriction enzyme used to cut the dna
        udna (str): 5' DNA string

    Returns:
        Iterator[Tuple[str, str]]: An iterator over the double stranded DNA segments
    """
    dnas: List[Iterator] = []
    for dna in ambiguous_dna_unambiguous(udna):
        seq = Seq(dna)
        fwd = map(str, enzyme.catalyse(seq))
        rev = map(
            lambda x: "".join(reversed(x)),
            reversed(enzyme.catalyse(seq.reverse_complement())),
        )

        dnas.append(iter(zip(fwd, rev)))
    return chain.from_iterable(dnas)


def find_overhang(
    fwd_len: int, rev_len: int, prev_fwd_hang: OverHang, prev_rev_hang: OverHang
) -> Tuple[OverHang, OverHang]:
    """
    Find DNA hanging ends.

    Takes care of four different cases:
        1. fwd_len == rev_len:
            a. ----------
               ----------
            b. ----------
                   ----------
        2. fwd_len != rev_len:
            a. ----------
               ------
            b. ----------
                   --
        And the reverse of these cases
    Args:
        fwd_len (int): Length of DNA string of 5' to 3'
        rev_len (int): Length of DNA string of 3' to 5'
        prev_fwd_hang (OverHang): previous 5' segment hanging lengths
        prev_rev_hang (OverHang): previous 3' segment hanging lengths

    Returns:
        Tuple[OverHang, OverHang]: A tuple with the lengths of the hanging DNA of both stands (5' and 3').
    """
    f_overhang = OverHang()
    t_overhang = OverHang()
    diff = fwd_len - rev_len
    if diff == 0:
        # In the case where both dna strands have the same length and the the hangs are equal
        # the new strands have no hangs so we can just use the default values of 0.
        if prev_fwd_hang.tprime != prev_rev_hang.fprime:
            if prev_fwd_hang.tprime > 0:
                f_overhang = f_overhang._replace(tprime=prev_fwd_hang.tprime)
            if prev_rev_hang.fprime > 0:
                f_overhang = f_overhang._replace(fprime=prev_rev_hang.fprime)
            t_overhang = t_overhang._replace(**f_overhang._asdict())
    else:
        if prev_fwd_hang.tprime == prev_rev_hang.fprime:
            if diff > 0:
                f_overhang = f_overhang._replace(tprime=diff)
            else:
                t_overhang = t_overhang._replace(fprime=-diff)
        else:
            f_overhang = f_overhang._replace(fprime=prev_rev_hang.fprime)
            t_overhang = t_overhang._replace(tprime=prev_fwd_hang.tprime)
            if diff > 0:
                if fwd_len - (rev_len + diff) > 0:
                    f_overhang = f_overhang._replace(tprime=fwd_len - (rev_len + diff))
            else:
                if rev_len - (fwd_len - diff) > 0:
                    t_overhang = t_overhang._replace(fprime=rev_len - (fwd_len - diff))
    return f_overhang, t_overhang


def find_overhand_dna(
    fwd_dna: str, rev_dna: str, fwd_overhang: OverHang, rev_overhang: OverHang
) -> Tuple[OverHangDNA, OverHangDNA]:
    """
    See ~find_overhang for cases and logic.

    Args:
        fwd_dna (str): The 5' DNA string
        rev_dna (str): The 3' DNA string:
        fwd_overhang (OverHang): The 5' strand overhangs
        rev_overhang (OverHang): The 5' strand overhangs

    Returns:
        Tuple[OverHangDNA, OverHangDNA]
    """
    return (
        OverHangDNA(
            fprime=fwd_dna[: fwd_overhang.fprime],
            tprime=fwd_dna[-fwd_overhang.tprime :] if fwd_overhang.tprime > 0 else "",
        ),
        OverHangDNA(
            tprime=rev_dna[: rev_overhang.tprime],
            fprime=rev_dna[-rev_overhang.fprime :] if rev_overhang.fprime > 0 else "",
        ),
    )
