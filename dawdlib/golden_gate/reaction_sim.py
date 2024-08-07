import collections
import json
from itertools import chain, combinations, permutations, product
from typing import (
    Any,
    Dict,
    Generator,
    Iterable,
    Iterator,
    List,
    NamedTuple,
    Optional,
    Set,
    Tuple,
    Union,
)

import networkx as nx
import numpy as np
import pandas as pd
from Bio.Restriction import Restriction
from Bio.Restriction.Restriction import Ov3, Ov5
from Bio.Seq import Seq
from tabulate import tabulate

from dawdlib.golden_gate.gate import Gate
from dawdlib.golden_gate.gate_data import GGData
from dawdlib.golden_gate.utils import (
    OligoTableEntry,
    Requirements,
    ambiguous_dna_unambiguous,
)


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
    is_const: bool = False

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


class ReactionGraphEdgeAttr(NamedTuple):
    """
    A structure for the edge attributes in reaction graph
    """

    id: str
    score: float


class ReactionGraph(nx.Graph):
    """
    A wrapper class around ~networkx.Graph which has a source and target fields.
    """

    source: DDNASection = DDNASection(
        fwd=SDNASection(start=-1, end=0), is_wt=True, is_const=True
    )
    target: DDNASection = DDNASection(
        fwd=SDNASection(start=np.iinfo(int).max, end=np.iinfo(int).max),
        is_wt=True,
        is_const=True,
    )


class ReactionGraphWt(NamedTuple):
    dna: str
    no_segments: int
    start: int
    end: int
    fidelity_sum: float


class ReactionSim:

    reaction_graph: ReactionGraph

    def __init__(self, ggdata: GGData, reqs: Requirements, enzymes: List[str]):
        self.ggdata = ggdata
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
                try:
                    dna_segs_it = enzyme_dna_segments(
                        enzyme, oligo_entry.full_oligo_dna
                    )
                except AssertionError as e:
                    e.args = (e.args, oligo_entry)
                    raise e
                start = (
                    oligo_entry.gate1.idx if not oligo_entry.gate1.src_or_target else 0
                )
                ddna_nodes.append(
                    get_ddna_sections(
                        enzyme,
                        dna_segs_it,
                        start,
                        oligo_entry.wt,
                        oligo_entry.const,
                    )
                )
        return chain.from_iterable(ddna_nodes)

    def create_reaction_graph(
        self, table_path: str
    ) -> Optional[Tuple[str, OligoTableEntry]]:
        """
        Creates the reaction graph given a table of ~dawdlib.goldengate.utils.OligoTableEntry
        Args:
            table_path (str): The path to a CSV file of OligoTableEntries.
        """
        try:
            ddna_iter = self.load_oligo_table(table_path)
        except AssertionError as e:
            return ", ".join(e.args[0]), OligoTableEntry(*e.args[1])
        rgraph = ReactionGraph()
        rgraph.add_nodes_from(ddna_iter)
        rgraph.add_edges_from(
            get_compatible_oligos(self.ggdata, rgraph.nodes, self.reqs.min_fidelity)
        )
        sources, targets = find_sources_targets(iter(rgraph.nodes))
        rgraph.add_edges_from(
            gen_src_trgt_edges(rgraph.source, rgraph.target, sources, targets)
        )
        self.reaction_graph = rgraph

    def get_wt_dna(self) -> Generator[Tuple[str, int, int, int, int], None, None]:
        """
        Finds all WT dna segments that can be connected in the reaction.

        Yields:
            str: The WT DNA sequence found in the reaction graph

        """
        sub_reaction_g = self.reaction_graph.subgraph(
            filter(lambda x: x.is_wt, self.reaction_graph.nodes)
        ).copy()
        pth: List[DDNASection]
        for pth in all_simple_paths(
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
            fidelity_sum = sum(
                map(lambda x: sub_reaction_g.edges[x]["score"], nx.utils.pairwise(pth))
            )
            yield ReactionGraphWt(dna, len(pth), start, end, fidelity_sum)

    def verify_reaction(
        self,
        expected_prod_len: int = 0,
        expected_no_segments: int = 0,
        expected_fidelity_sum: float = 0,
    ) -> Tuple[bool, int, List[DDNASection]]:
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
        sub_reaction_g = self.reaction_graph.subgraph(
            filter(lambda x: x.is_wt == x.is_const, self.reaction_graph.nodes)
        ).copy()
        pth: List[DDNASection]
        cnt = -1
        for cnt, pth in enumerate(
            all_simple_paths(
                sub_reaction_g,
                self.reaction_graph.source,
                self.reaction_graph.target,
            )
        ):
            if not all(
                (
                    pth[i].fwd.start < pth[i + 1].fwd.start
                    and pth[i].fwd.end < pth[i + 1].fwd.end
                    for i in range(1, len(pth) - 1)
                )
            ):
                return False, cnt + 1, pth
            if 0 < expected_prod_len != len("".join(node.fwd.dna for node in pth)):
                return False, cnt + 1, pth
            if 0 < expected_no_segments != len(pth):
                return False, cnt + 1, pth
            if (
                0
                < expected_fidelity_sum
                != sum(
                    map(
                        lambda x: sub_reaction_g.edges[x]["score"],
                        nx.utils.pairwise(pth),
                    )
                )
            ):
                return False, cnt + 1, pth
        return True, cnt + 1, []

    def get_all_products(self) -> Generator[Tuple[str, int, int, int], None, None]:
        """
        Use wisely! this method yields all possible products of the golden gate reaction.
        Yields:
            Tuple: [dna, number of segments, start position, end position]
        """
        pth: List[DDNASection]
        for pth in all_simple_paths(
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
        min_fidelity: float,
        neb_table_temp: int = 25,
        neb_table_time: int = 18,
    ):

        self.sim = ReactionSim(
            GGData(temperature=neb_table_temp, hours=neb_table_time),
            Requirements(
                0, 0, 0, min_fidelity=min_fidelity, min_efficiency=GGData.MIN_EFFICIENCY
            ),
            enzymes,
        )
        self.sim.create_reaction_graph(table_path)
        return self

    def get_wt(
        self,
        table_path: str,
        enzymes: Union[List[str], str],
        min_fidelity: float,
        neb_table_temp: int = 25,
        neb_table_time: int = 18,
    ):
        """
        Prints the WT sequence(s) created by "simulating" the golden gate reaction.

        Args:

            table_path (str): Path to segments table
            enzymes (str|list[str]): which enzyme(s) to use in the reaction given. Enzymes should be given in a list format, i.e. [BsaI]. You must respect the usual naming convention with the upper case letters and Latin numbering (in upper case as well).
            min_fidelity (float) : The fidelity above overhangs are considered connected
            neb_table_temp (int): Specifies which NEB data to use (default=25)
            neb_table_time (int): Specifies which NEB data to use (default=18)
        """
        if isinstance(enzymes, str):
            enzymes = [enzymes]
        self._setup(table_path, enzymes, min_fidelity, neb_table_temp, neb_table_time)
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
        min_fidelity: float,
        expected_dna_len: int = 0,
        expected_no_segments: int = 0,
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
            min_fidelity (float) : The fidelity above overhangs are considered connected
            expected_dna_len (int): The number of BPs each product is expected to contain
            expected_no_segments (int): The number of segments each product is expected to contain
            neb_table_temp (int): Specifies which NEB data to use (default=25)
            neb_table_time (int): Specifies which NEB data to use (default=18)

        """
        if isinstance(enzymes, str):
            enzymes = [enzymes]
        self._setup(table_path, enzymes, min_fidelity, neb_table_temp, neb_table_time)
        success, count, extra = self.sim.verify_reaction(
            expected_prod_len=expected_dna_len,
            expected_no_segments=expected_no_segments,
        )
        if success:
            print(
                f"Reaction simulation is successful, {count} unique sequences were found."
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
        min_fidelity: float,
        neb_table_temp: int = 25,
        neb_table_time: int = 18,
    ):
        """
        Use with caution!!!
        This method prints all golden gate reaction products in a table format.
        Args:
            table_path (str): Path to segments table
            enzymes (str|list[str]): which enzyme(s) to use in the reaction given. Enzymes should be given in a list format, i.e. [BsaI]. You must respect the usual naming convention with the upper case letters and Latin numbering (in upper case as well).
            min_fidelity (float) : The fidelity above overhangs are considered connected
            neb_table_temp (int): Specifies which NEB data to use (default=25)
            neb_table_time (int): Specifies which NEB data to use (default=18)
        """
        if isinstance(enzymes, str):
            enzymes = [enzymes]
        self._setup(table_path, enzymes, min_fidelity, neb_table_temp, neb_table_time)
        products = self.sim.get_all_products()
        print(tabulate(products, headers=["DNA", "No. segments", "Start", "End"]))


def gen_src_trgt_edges(
    source: DDNASection,
    target: DDNASection,
    sources: Iterable[DDNASection],
    targets: Iterable[DDNASection],
) -> Iterator[Tuple[DDNASection, DDNASection, Dict[str, Any]]]:
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
    src_atrr = ReactionGraphEdgeAttr(id="SRC", score=0)
    target_atrr = ReactionGraphEdgeAttr(id="TARGET", score=0)
    return chain(
        map(lambda x: x + (src_atrr._asdict(),), product([source], sources)),
        map(lambda x: x + (target_atrr._asdict(),), product([target], targets)),
    )


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
    ggdata: GGData, ddnas: Iterable[DDNASection], min_fidelity: float
) -> Generator[Tuple[DDNASection, DDNASection, Dict], None, None]:
    """

    Args:
        ggdata (GGData): A wrapper around the gates data
        ddnas (Iterator[DDNASection]):
        min_fidelity (float): fidelity above this score marks the two dna sections as compatible

    Yields:
         Tuple[DDNASection, DDNASection] - Signifying these two dna sections can connect to one another in the reaction
    """
    d1: DDNASection
    d2: DDNASection
    for d1, d2 in combinations(ddnas, 2):
        fidelity = max(
            chain.from_iterable(
                map(
                    lambda x: (
                        ggdata.fwd_fidelity.at[x],
                        ggdata.fwd_fidelity.at[x[::-1]],
                    ),
                    product(
                        filter(None, d1.get_hang_dna_rev()),
                        filter(None, d2.get_hang_dna_rev()),
                    ),
                )
            )
        )
        if fidelity >= 10:
            rgea = ReactionGraphEdgeAttr(
                id="-".join(
                    map(
                        lambda x: f"({x[0]},{x[1]})",
                        sorted(
                            [(d1.fwd.start, d1.fwd.end), (d2.fwd.start, d2.fwd.end)]
                        ),
                    )
                ),
                score=fidelity,
            )
            yield d1, d2, rgea._asdict()


def get_ddna_sections(
    enzyme: Union[Ov5, Ov3],
    dna_segs: Iterator[Tuple[str, str]],
    start: int,
    is_wt: bool = False,
    is_const: bool = False,
) -> Generator[DDNASection, None, None]:
    """

    Args:
        enzyme (Ov5): Restriction enzyme which was used to digest the DNA.
        dna_segs (Iterator[Tuple[str, str]]): An iterator over double stranded DNA segments
                                              (after digestion by restrction enzyme).
        start (int): An index pointing to the start of the DNA segment or the location of the gate.
        is_wt (bool): Whether the current oligo is WT or not.
        is_const (bool): Whether the current oligo is constant or not.

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
            is_const=is_const,
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
        assert (
            len(enzyme.search(seq)) < 3
        ), "Dna has more than two occurrences of restriction site!!!"
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


def all_simple_paths(
    G: nx.Graph, source: DDNASection, target: DDNASection, cutoff: Optional[int] = None
) -> Iterable[List[DDNASection]]:
    """Generate all simple paths in the graph G from source to target.

    A simple path is a path with no repeated nodes.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path

    target : nodes
       Single node or iterable of nodes at which to end path

    cutoff : integer, optional
        Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    path_generator: generator
       A generator that produces lists of simple paths.  If there are no paths
       between the source and target within the given cutoff the generator
       produces no output.

    Examples
    --------
    This iterator generates lists of nodes::

        >>> G = nx.complete_graph(4)
        >>> for path in nx.all_simple_paths(G, source=0, target=3):
        ...     print(path)
        ...
        [0, 1, 2, 3]
        [0, 1, 3]
        [0, 2, 1, 3]
        [0, 2, 3]
        [0, 3]

    You can generate only those paths that are shorter than a certain
    length by using the `cutoff` keyword argument::

        >>> paths = nx.all_simple_paths(G, source=0, target=3, cutoff=2)
        >>> print(list(paths))
        [[0, 1, 3], [0, 2, 3], [0, 3]]

    To get each path as the corresponding list of edges, you can use the
    :func:`networkx.utils.pairwise` helper function::

        >>> paths = nx.all_simple_paths(G, source=0, target=3)
        >>> for path in map(nx.utils.pairwise, paths):
        ...     print(list(path))
        [(0, 1), (1, 2), (2, 3)]
        [(0, 1), (1, 3)]
        [(0, 2), (2, 1), (1, 3)]
        [(0, 2), (2, 3)]
        [(0, 3)]

    Pass an iterable of nodes as target to generate all paths ending in any of several nodes::

        >>> G = nx.complete_graph(4)
        >>> for path in nx.all_simple_paths(G, source=0, target=[3, 2]):
        ...     print(path)
        ...
        [0, 1, 2]
        [0, 1, 2, 3]
        [0, 1, 3]
        [0, 1, 3, 2]
        [0, 2]
        [0, 2, 1, 3]
        [0, 2, 3]
        [0, 3]
        [0, 3, 1, 2]
        [0, 3, 2]

    Iterate over each path from the root nodes to the leaf nodes in a
    directed acyclic graph using a functional programming approach::

        >>> from itertools import chain
        >>> from itertools import product
        >>> from itertools import starmap
        >>> from functools import partial
        >>>
        >>> chaini = chain.from_iterable
        >>>
        >>> G = nx.DiGraph([(0, 1), (1, 2), (0, 3), (3, 2)])
        >>> roots = (v for v, d in G.in_degree() if d == 0)
        >>> leaves = (v for v, d in G.out_degree() if d == 0)
        >>> all_paths = partial(nx.all_simple_paths, G)
        >>> list(chaini(starmap(all_paths, product(roots, leaves))))
        [[0, 1, 2], [0, 3, 2]]

    The same list computed using an iterative approach::

        >>> G = nx.DiGraph([(0, 1), (1, 2), (0, 3), (3, 2)])
        >>> roots = (v for v, d in G.in_degree() if d == 0)
        >>> leaves = (v for v, d in G.out_degree() if d == 0)
        >>> all_paths = []
        >>> for root in roots:
        ...     for leaf in leaves:
        ...         paths = nx.all_simple_paths(G, root, leaf)
        ...         all_paths.extend(paths)
        >>> all_paths
        [[0, 1, 2], [0, 3, 2]]

    Iterate over each path from the root nodes to the leaf nodes in a
    directed acyclic graph passing all leaves together to avoid unnecessary
    compute::

        >>> G = nx.DiGraph([(0, 1), (2, 1), (1, 3), (1, 4)])
        >>> roots = (v for v, d in G.in_degree() if d == 0)
        >>> leaves = [v for v, d in G.out_degree() if d == 0]
        >>> all_paths = []
        >>> for root in roots:
        ...     paths = nx.all_simple_paths(G, root, leaves)
        ...     all_paths.extend(paths)
        >>> all_paths
        [[0, 1, 3], [0, 1, 4], [2, 1, 3], [2, 1, 4]]

    Notes
    -----
    This algorithm uses a modified depth-first search to generate the
    paths [1]_.  A single path can be found in $O(V+E)$ time but the
    number of simple paths in a graph can be very large, e.g. $O(n!)$ in
    the complete graph of order $n$.

    References
    ----------
    .. [1] R. Sedgewick, "Algorithms in C, Part 5: Graph Algorithms",
       Addison Wesley Professional, 3rd ed., 2001.

    """
    if source not in G:
        raise nx.NodeNotFound("source node %s not in graph" % str(source))
    if target not in G:
        raise nx.NodeNotFound("target node %s not in graph" % str(target))
    targets = {target}
    if source in targets:
        return []
    if cutoff is None:
        cutoff = len(G) - 1
    if cutoff < 1:
        return []
    if G.is_multigraph():
        raise ValueError("Not implemented for multi-graphs")
    else:
        return _all_simple_paths_graph(G, source, targets, cutoff)


def _all_simple_paths_graph(
    G: nx.Graph, source: DDNASection, targets: Set[DDNASection], cutoff: int
) -> Generator[List[DDNASection], None, None]:
    class Child(NamedTuple):
        incoming_id: str
        child: DDNASection

    visited: Dict[DDNASection, Any] = collections.OrderedDict.fromkeys([source])
    stack: List[Iterator[Child]] = [
        map(lambda x: Child(incoming_id=x[1]["id"], child=x[0]), G[source].items())
    ]
    while stack:
        children = stack[-1]
        incoming_id, child = next(children, (None, None))
        if child is None:
            stack.pop()
            visited.popitem()
        elif len(visited) < cutoff:
            if child in visited:
                continue
            if child in targets:
                yield list(visited) + [child]
            visited[child] = None
            if targets - set(visited.keys()):  # expand stack until find all targets
                stack.append(
                    filter(
                        lambda x: x.incoming_id != incoming_id,
                        map(
                            lambda x: Child(incoming_id=x[1]["id"], child=x[0]),
                            G[child].items(),
                        ),
                    )
                )
            else:
                visited.popitem()  # maybe other ways to child
        else:  # len(visited) == cutoff:
            for target in (targets & (set(children) | {child})) - set(visited.keys()):
                yield list(visited) + [target]
            stack.pop()
            visited.popitem()
