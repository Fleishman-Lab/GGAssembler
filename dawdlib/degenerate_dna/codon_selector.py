import typing as tp
from itertools import chain, combinations, product

import networkx as nx
import numpy as np

from dawdlib.degenerate_dna.alg_x import AlgXSolver
from dawdlib.degenerate_dna.codon_utils import AmbigCodon
from dawdlib.degenerate_dna.codon_utils import CodonSelector as UtilsCodonSelector
from dawdlib.degenerate_dna.codon_utils import _get_score
from dawdlib.degenerate_dna.structs import DEG_NUCL_CODES, PosCodon


class CodonSelector:
    """Usage:
    initiate giving an organism id as an NCBI Taxonomy id.
    Commons ids:
        Escherichia coli ID = "37762"
        Saccharomyces cerevisiae ID = "4932"

    Lookup ids at: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?

    For paper see: https://doi.org/10.1093/nar/gkr1178

    To get best degenerate codon library call optimise_codons with the
    required amino acids as single letter list.
    example optimise_codons(['L', 'D'])
    """

    def __init__(self, organism_id: str):
        self.organism_id = organism_id
        self._graph: tp.Optional[nx.Graph] = None
        self._codon_dict: tp.Optional[tp.Dict] = None
        self.cod_sel = UtilsCodonSelector()

    def _set_graph(self, graph: nx.Graph):
        self._graph = graph

    def _get_graph(self) -> tp.Optional[nx.Graph]:
        return self._graph

    graph = property(_get_graph, _set_graph)

    def _set_codon_dict(self, codon_dict: tp.Dict[str, tp.Dict]):
        self._codon_dict = codon_dict

    def _get_codon_dict(self) -> tp.Optional[tp.Dict]:
        return self._codon_dict

    codon_dict = property(_get_codon_dict, _set_codon_dict)

    def select_aminoacids_codon(self, amino_acids: str) -> tp.Iterable[AmbigCodon]:
        req_aas = set(list(amino_acids))

        for codon in self.cod_sel.optimise_codons(amino_acids, self.organism_id):
            codon_aas = {amino_acid.amino_acid for amino_acid in codon.amino_acids}
            if req_aas.issuperset(codon_aas):
                yield codon

    def select_codon_aminoacids(self, amino_acids: str) -> tp.Iterable[AmbigCodon]:
        req_aas = set(list(amino_acids))
        for bps in product(DEG_NUCL_CODES, repeat=3):
            for codon in self.cod_sel.analyse_codon("".join(bps), self.organism_id):
                codon_aas = {amino_acid.amino_acid for amino_acid in codon.amino_acids}
                if req_aas.issuperset(codon_aas):
                    codon._replace(
                        amino_acids=tuple(
                            map(lambda x: x._replace(type=1), codon.amino_acids)
                        )
                    )
                    yield codon

    def _gen_degenerate_codons_big(
        self, amino_acids: tp.List[str]
    ) -> tp.List[PosCodon]:
        res = self.select_codon_aminoacids("".join(amino_acids))
        res2 = _codon_to_codons_gen(res)
        return list(res2)

    def _gen_degenerate_codons_sm(self, amino_acids: tp.List[str]) -> tp.List[PosCodon]:
        res = (
            self.select_aminoacids_codon("".join(amino_acid_set))
            for amino_acid_set in _powerset(amino_acids)
        )
        res2 = (_codon_to_codons_gen(x) for x in res)
        return list(chain(*res2))

    def _gen_degenerate_codons(self, amino_acids: tp.List[str]) -> tp.List[PosCodon]:
        if 2 ** len(amino_acids) > 3375:
            return self._gen_degenerate_codons_big(amino_acids)
        return self._gen_degenerate_codons_sm(amino_acids)

    def _lookup_degenerate_codons(self, amino_acids: tp.List[str]) -> tp.List[PosCodon]:
        acc = []
        for amino_acid_set in _powerset(amino_acids):
            try:
                acc.append(self.codon_dict["".join(sorted(amino_acid_set))])
            except KeyError:
                pass
        return acc

    def optimise_codons(self, amino_acids: tp.List[str], mode="exact") -> PosCodon:

        if mode == "graph" and self.graph is not None:
            return _optimise_codons_graph(amino_acids, graph=self.graph)

        if self.codon_dict is not None:
            codons = self._lookup_degenerate_codons(amino_acids)
        else:
            codons = self._gen_degenerate_codons(amino_acids)
            codons = _remove_duplicate_codons(codons)
        codons = sorted(codons, key=lambda x: len(x.encoded_acids), reverse=True)

        if mode == "graph":
            return _optimise_codons_graph(amino_acids, codons)
        if mode == "greedy":
            return _optimise_codons_greedy(amino_acids, codons)
        if mode == "exact":
            return _optimise_codons_exact(amino_acids, codons)

        raise ValueError(
            "Unknown value %s for mode given.\
                Allowed value are: lp, graph, greedy, exact"
            % mode
        )


P = tp.TypeVar("P")


def _powerset(p: tp.Sequence[P]) -> tp.Generator[tp.List[P], None, None]:
    for i in range(1, len(p) + 1):
        for item in combinations(p, i):
            yield list(item)


def _codon_to_codons_gen(codon_gen: tp.Iterable[AmbigCodon]) -> tp.Iterable[PosCodon]:
    for codon in codon_gen:
        yield _codon_to_codons(codon)


def _codon_to_codons(codon: AmbigCodon) -> PosCodon:
    score = codon.score
    if 0 == score:
        score = _get_score(codon.amino_acids)

    return PosCodon(
        ambiguous_codons=(codon.ambiguous_codon,),
        amino_acids=codon.amino_acids,
        encoded_acids=frozenset(
            amino_acid.amino_acid for amino_acid in codon.amino_acids
        ),
        score=score,
    )


def _combine_condons(codons_list: tp.Sequence[PosCodon]) -> PosCodon:
    cdn = PosCodon(
        ambiguous_codons=tuple(
            chain(*(codons.ambiguous_codons for codons in codons_list))
        ),
        amino_acids=tuple(chain(*(codons.amino_acids for codons in codons_list))),
        encoded_acids=frozenset(
            (aa for codons in codons_list for aa in codons.encoded_acids)
        ),
    )
    return cdn._replace(score=_get_score(cdn.amino_acids))


def _remove_duplicate_codons(codons_list: tp.Sequence[PosCodon]) -> tp.List[PosCodon]:
    codon_dict: tp.Dict[str, PosCodon] = {}
    for codon in codons_list:
        try:
            if codon_dict["".join(sorted(codon.encoded_acids))].score < codon.score:
                codon_dict["".join(sorted(codon.encoded_acids))] = codon
        except KeyError:
            codon_dict["".join(sorted(codon.encoded_acids))] = codon
    return list(codon_dict.values())


def _is_codon_set_valid(
    req_aas: set, codon_aas: tp.Union[tp.Set, tp.FrozenSet]
) -> bool:
    return len(req_aas.symmetric_difference(codon_aas)) == 0


def _add_graph_edges(
    g: nx.Graph, nodes: tp.List[str], edges: tp.List[str], codons: tp.List[tp.Dict]
):
    for v1, v2 in combinations(nodes, 2):
        diff = "".join(sorted(set(v1) - set(v2)))
        try:
            idx = edges.index(diff)
            g.add_edge(v2, v1, codon=codons[idx])
        except ValueError:
            pass


def _optimise_codons_graph(
    amino_acids: tp.List[str],
    codons: tp.Optional[tp.List[PosCodon]] = None,
    graph: tp.Optional[nx.Graph] = None,
) -> PosCodon:
    """Build a graph of possible degenerate codons
        that only encode the required given codons.

        The nodes are all possible combinations of amino acids
        given. (No. of vertices is 2 **|amino acids|)

        Edge are degenerate codons which encode the required AAs.
        Two nodes are connected by an edge if there is a degenerate codon
        that can create the difference between the two nodes AAs.
        """
    best_codon: PosCodon = PosCodon()
    best_codon_score = np.NINF
    aa_set = frozenset(amino_acids)

    if codons is not None:
        g = _generate_codons_set_graph(codons)
    else:
        g = graph

    for path in nx.all_shortest_paths(g, frozenset(), aa_set):
        codons_list = [
            PosCodon(**g[v1][v2]["codon"]) for v1, v2 in zip(path[:-1], path[1:])
        ]
        codon_comb = _combine_condons(codons_list)
        if codon_comb.score > best_codon_score:
            best_codon = codon_comb
            best_codon_score = codon_comb.score

    return best_codon


def _generate_codons_graph(
    amino_acids: tp.List[str], codons: tp.List[tp.Dict]
) -> nx.Graph:
    nodes = ["".join(sorted(x)) for x in _powerset(amino_acids)]
    nodes.append("")
    nodes = sorted(nodes, key=len, reverse=True)
    edges = ["".join(sorted(x["encoded_acids"])) for x in codons]
    g = nx.Graph()
    g.add_nodes_from(nodes)
    _add_graph_edges(g, nodes, edges, codons)
    return g


def _generate_codons_set_graph(codons: tp.List[PosCodon]) -> nx.Graph:

    stack: tp.List[frozenset] = [frozenset()]

    nodes: tp.List[tp.Tuple[frozenset, tp.Dict]] = [(frozenset(), PosCodon()._asdict())]
    edges: tp.List[tp.Tuple[frozenset, frozenset, tp.Dict]] = []

    visited = set()

    while stack:
        u = stack.pop()
        visited.add(u)
        for codon in codons:
            encoded_acids = codon.encoded_acids
            if u.isdisjoint(encoded_acids):
                v = u.union(encoded_acids)
                edges.append((u, v, {"codon": codon._asdict()}))
                if v not in visited:
                    nodes.append((v, codon._asdict()))
                    stack.append(v)

    g = nx.Graph()
    g.add_nodes_from(nodes)
    g.add_edges_from(edges)
    return g


def _generate_reverse_codons_graph(
    amino_acids: tp.List[str], codons: tp.List[tp.Dict]
) -> tp.Tuple[nx.Graph, tp.Tuple[frozenset, dict], tp.Tuple[frozenset, dict]]:
    nodes = [(codon["encoded_acids"], codon) for codon in codons]
    source: tp.Tuple[frozenset, tp.Dict] = (frozenset(), {})
    target: tp.Tuple[frozenset, tp.Dict] = (frozenset(amino_acids), {})
    nodes += [source, target]
    edges = [(v1, v2) for v1, v2 in combinations(nodes, 2) if len(v1[0] & v2[0]) == 0]

    g = nx.Graph()
    g.add_nodes_from(nodes)
    g.add_edges_from(edges)
    return g, source, target


def _optimise_codons_greedy(
    amino_acids: tp.List[str], codons: tp.List[PosCodon]
) -> PosCodon:

    req_aas = set(amino_acids)
    not_found_aas = set(amino_acids)
    used_codons: tp.List[PosCodon] = []
    codons = sorted(codons, key=lambda x: len(x.encoded_acids), reverse=False)

    while not_found_aas and codons:
        cdn = codons.pop()
        used_codons.append(cdn)
        not_found_aas = not_found_aas.difference(cdn.encoded_acids)
        codons = filter(lambda x: not_found_aas.intersection(x.encoded_acids), codons)
        codons = sorted(
            codons,
            key=lambda x: (
                len(not_found_aas.intersection(x.encoded_acids)),
                -len(x.encoded_acids.difference(not_found_aas)),
            ),
            reverse=False,
        )

    codon_comb = _combine_condons(used_codons)
    if _is_codon_set_valid(req_aas, codon_comb.encoded_acids):
        return codon_comb

    raise ValueError("No codons were found for AAs %s" % amino_acids)


def _optimise_codons_exact(
    amino_acids: tp.List[str], codons: tp.List[PosCodon]
) -> PosCodon:
    best_codon: PosCodon = PosCodon()
    best_codon_score = np.NINF
    best_codon_len = np.inf
    req_aas = set(amino_acids)

    for cdns in AlgXSolver.solve(req_aas, codons):
        if len(cdns) > best_codon_len:
            continue
        codon_comb = _combine_condons(cdns)
        if not _is_codon_set_valid(req_aas, codon_comb.encoded_acids):
            continue
        if len(codon_comb.ambiguous_codons) < best_codon_len:
            best_codon = codon_comb
            best_codon_len = len(codon_comb.ambiguous_codons)
            best_codon_score = codon_comb.score
            continue
        # Last case: lengths are equal
        if codon_comb.score > best_codon_score:
            best_codon = codon_comb
            best_codon_len = len(codon_comb.ambiguous_codons)
            best_codon_score = codon_comb.score

    if _is_codon_set_valid(req_aas, best_codon.encoded_acids):
        return best_codon

    raise ValueError("No codons were found for AAs %s" % amino_acids)


def _populate_set_matrix(
    a: np.ndarray, amino_acids: tp.List[str], codons: tp.List[tp.Dict]
):
    for i, codon in enumerate(codons):
        for aa in codon["encoded_acids"]:
            a[amino_acids.index(aa), i] = 1
