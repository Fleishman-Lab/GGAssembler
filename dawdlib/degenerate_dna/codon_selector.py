import typing as tp
from itertools import chain, combinations, product

import cvxpy as cp
import networkx as nx
import numpy as np
from dawdlib.degenerate_dna.codon_utils import CodonSelector as UtilsCodonSelector
from dawdlib.degenerate_dna.codon_utils import _get_score

DEG_NUCL_CODES = [
    "A",
    "C",
    "G",
    "T",
    "R",
    "Y",
    "S",
    "W",
    "K",
    "M",
    "B",
    "D",
    "H",
    "V",
    "N",
]


class CodonSelector:
    """Usage:
    initiate giving an organism id as an NCBI Taxonomy id
        See: https://doi.org/10.1093/nar/gkr1178

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

    def select_aminoacids_codon(
        self, amino_acids: str
    ) -> tp.Generator[tp.Any, None, None]:
        for codon in self.cod_sel.optimise_codons(amino_acids, self.organism_id):
            codon_aas = {
                amino_acid["amino_acid"] for amino_acid in codon["amino_acids"]
            }
            if set(list(amino_acids)).issuperset(set(codon_aas)):
                yield codon

    def select_codon_aminoacids(
        self, amino_acids: str
    ) -> tp.Generator[tp.Any, None, None]:
        for bps in product(DEG_NUCL_CODES, repeat=3):
            for codon in self.cod_sel.analyse_codon("".join(bps), self.organism_id):
                codon_aas = {
                    amino_acid["amino_acid"] for amino_acid in codon["amino_acids"]
                }
                if set(list(amino_acids)).issuperset(set(codon_aas)):
                    for amino_acid in codon["amino_acids"]:
                        amino_acid["type"] = 1
                    yield codon

    def _gen_degenerate_codons_big(self, amino_acids: tp.List[str]) -> tp.List[tp.Dict]:
        res = self.select_codon_aminoacids("".join(amino_acids))
        res2 = _codon_to_codons_gen(res)
        return list(res2)

    def _gen_degenerate_codons_sm(self, amino_acids: tp.List[str]) -> tp.List[tp.Dict]:
        res = (
            self.select_aminoacids_codon("".join(amino_acid_set))
            for amino_acid_set in _powerset(amino_acids)
        )
        res2 = (_codon_to_codons_gen(x) for x in res)
        return list(chain(*res2))

    def _gen_degenerate_codons(self, amino_acids: tp.List[str]) -> tp.List[tp.Dict]:
        if 2 ** len(amino_acids) > 3375:
            return self._gen_degenerate_codons_big(amino_acids)
        return self._gen_degenerate_codons_sm(amino_acids)

    def _lookup_degenerate_codons(self, amino_acids: tp.List[str]) -> tp.List[tp.Dict]:
        acc = []
        for amino_acid_set in _powerset(amino_acids):
            try:
                acc.append(self.codon_dict["".join(sorted(amino_acid_set))])
            except KeyError:
                pass
        return acc

    def optimise_codons(self, amino_acids: tp.List[str], mode="graph") -> tp.Dict:

        if mode == "graph" and self.graph is not None:
            return _optimise_codons_graph(amino_acids, graph=self.graph)

        if self.codon_dict is not None:
            codons = self._lookup_degenerate_codons(amino_acids)
        else:
            codons = self._gen_degenerate_codons(amino_acids)
            codons = _remove_duplicate_codons(codons)
        codons = sorted(codons, key=lambda x: len(x["encoded_acids"]), reverse=True)

        if mode == "lp":
            return _optimise_codons_lp(amino_acids, codons)
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


def _codon_to_codons_gen(
    codon_gen: tp.Generator[tp.Any, None, None]
) -> tp.Generator[tp.Dict, None, None]:
    for codon in codon_gen:
        yield _codon_to_codons(codon)


def _codon_to_codons(codon: tp.Dict) -> tp.Dict:
    try:
        score = codon["score"]
    except KeyError:
        score = _get_score(codon["amino_acids"])

    return {
        "ambiguous_codons": [codon["ambiguous_codon"]],
        "amino_acids": codon["amino_acids"],
        "encoded_acids": frozenset(
            amino_acid["amino_acid"] for amino_acid in codon["amino_acids"]
        ),
        "score": score,
    }


def _combine_condons(codons_list: tp.Sequence[tp.Dict]) -> tp.Dict:
    codons_dict = {
        "ambiguous_codons": list(
            chain(*(codons["ambiguous_codons"] for codons in codons_list))
        ),
        "amino_acids": list(chain(*(codons["amino_acids"] for codons in codons_list))),
        "encoded_acids": frozenset(
            (aa for codons in codons_list for aa in codons["encoded_acids"])
        ),
    }
    codons_dict["score"] = _get_score(codons_dict["amino_acids"])
    return codons_dict


def _remove_duplicate_codons(codons_list: tp.Sequence[tp.Dict]) -> tp.List[tp.Dict]:
    codon_dict: tp.Dict[str, tp.Dict] = {}
    for codon in codons_list:
        try:
            if (
                codon_dict["".join(sorted(codon["encoded_acids"]))]["score"]
                < codon["score"]
            ):
                codon_dict["".join(sorted(codon["encoded_acids"]))] = codon
        except KeyError:
            codon_dict["".join(sorted(codon["encoded_acids"]))] = codon
    return list(codon_dict.values())


def _is_codon_set_valid(req_aas: set, codon_aas: set) -> bool:
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
    codons: tp.Optional[tp.List[tp.Dict]] = None,
    graph: tp.Optional[nx.Graph] = None,
) -> tp.Dict:
    """Build a graph of possible degenerate codons
        that only encode the required given codons.

        The nodes are all possible combinations of amino acids
        given. (No. of vertices is 2 **|amino acids|)

        Edge are degenerate codons which encode the required AAs.
        Two nodes are connected by an edge if there is a degenerate codon
        that can create the difference between the two nodes AAs.
        """
    best_codon: tp.Dict = {}
    best_codon_score = np.NINF
    # aa_str = "".join(sorted(amino_acids))
    aa_set = frozenset(amino_acids)

    if codons is not None:
        g = _generate_codons_set_graph(codons)
    else:
        g = graph

    for path in nx.all_shortest_paths(g, frozenset(), aa_set):
        codons_list = [g[v1][v2]["codon"] for v1, v2 in zip(path[:-1], path[1:])]
        codon_comb = _combine_condons(codons_list)
        if codon_comb["score"] > best_codon_score:
            best_codon = codon_comb
            best_codon_score = codon_comb["score"]

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


def _generate_codons_set_graph(codons: tp.List[tp.Dict]) -> nx.Graph:

    stack: tp.List[frozenset] = [frozenset()]

    nodes: tp.List[tp.Tuple[frozenset, tp.Dict]] = [(frozenset(), {})]
    edges: tp.List[tp.Tuple[frozenset, frozenset, tp.Dict]] = []

    visited = set()

    while stack:
        u = stack.pop()
        visited.add(u)
        for codon in codons:
            encoded_acids = codon["encoded_acids"]
            if u.isdisjoint(encoded_acids):
                v = u.union(encoded_acids)
                edges.append((u, v, {"codon": codon}))
                if v not in visited:
                    nodes.append((v, codon))
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
    amino_acids: tp.List[str], codons: tp.List[tp.Dict]
) -> tp.Dict:
    req_aas_set = set(amino_acids)

    for codons_list in _powerset(codons):
        codon_comb = _combine_condons(codons_list)
        if _is_codon_set_valid(req_aas_set, codon_comb["encoded_acids"]):
            return codon_comb
    return {}


def _optimise_codons_exact(
    amino_acids: tp.List[str], codons: tp.List[tp.Dict]
) -> tp.Dict:
    best_codon: tp.Dict = {}
    best_codon_score = np.NINF
    best_codon_len = np.inf
    req_aas_set = set(amino_acids)

    for codons_list in _powerset(codons):
        codon_comb = _combine_condons(codons_list)
        if not _is_codon_set_valid(req_aas_set, codon_comb["encoded_acids"]):
            continue
        if len(codon_comb["ambiguous_codons"]) > best_codon_len:
            continue
        if codon_comb["score"] < best_codon_score:
            continue
        best_codon = codon_comb
        best_codon_len = len(codon_comb["ambiguous_codons"])
        best_codon_score = codon_comb["score"]
    return best_codon


def _populate_set_matrix(
    a: np.array, amino_acids: tp.List[str], codons: tp.List[tp.Dict]
) -> None:
    for i, codon in enumerate(codons):
        for aa in codon["encoded_acids"]:
            a[amino_acids.index(aa), i] = 1


def _optimise_codons_lp(amino_acids: tp.List[str], codons: tp.List[tp.Dict]) -> tp.Dict:

    x = cp.Variable(len(codons), boolean=True)
    a = np.zeros((len(amino_acids), len(codons)))
    _populate_set_matrix(a, amino_acids, codons)
    req = np.ones((len(amino_acids),))
    constraints = [a @ x >= req]
    objective = cp.Minimize(cp.sum(x))
    problem = cp.Problem(objective, constraints=constraints)
    problem.solve()
    selected_codons = [codons[i] for i, val in enumerate(x.value) if val >= 0.5]
    return _combine_condons(selected_codons)
