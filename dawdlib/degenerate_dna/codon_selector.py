import typing as tp
from itertools import chain, combinations

import networkx as nx
import numpy as np
from codon_utils import CodonSelector as UtilsCodonSelector
from codon_utils import _get_score


class CodonSelector(object):
    '''Usage:
    initiate giving an organism id as an NCBI Taxonomy id
        See: https://doi.org/10.1093/nar/gkr1178

    To get best degenerate codon library call optimise_codons with the
    required amino acids as single letter list.
    example optimise_codons(['L', 'D'])
    '''
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

    def select_codons(self,
                      amino_acids: str) -> tp.Generator[tp.Any, None, None]:
        for codon in self.cod_sel.optimise_codons(amino_acids,
                                                  self.organism_id):
            codon_aas = set([
                amino_acid['amino_acid'] for amino_acid in codon['amino_acids']
            ])
            if set(list(amino_acids)).issuperset(set(codon_aas)):
                yield codon

    def _gen_degenerate_codons(self,
                               amino_acids: tp.List[str]) -> tp.List[tp.Dict]:
        res = (self.select_codons(''.join(amino_acid_set))
               for amino_acid_set in _powerset(amino_acids))
        res2 = (_codon_to_codons_gen(x) for x in res)
        return list(chain(*res2))

    def _lookup_degenerate_codons(self, amino_acids: tp.List[str]
                                  ) -> tp.List[tp.Dict]:
        acc = []
        for amino_acid_set in _powerset(amino_acids):
            try:
                acc.append(self.codon_dict[''.join(sorted(amino_acid_set))])
            except KeyError:
                pass
        return acc

    def _optimise_codons_graph(self, amino_acids: tp.List[str],
                               codons: tp.Optional[tp.List[tp.Dict]]
                               ) -> tp.Optional[tp.Dict]:
        '''Build a graph of possible degenerate codons
        that only encode the required given codons.

        The nodes are all possible combinations of amino acids
        given. (No. of vertices is 2 **|amino acids|)

        Edge are degenerate codons which encode the required AAs.
        Two nodes are connected by an edge if there is a degenerate codon
        that can create the difference between the two nodes AAs.
        '''
        best_codon = None
        best_codon_score = np.NINF
        aa_str = ''.join(sorted(amino_acids))

        if codons is not None:
            nodes = [''.join(sorted(x)) for x in _powerset(amino_acids)]
            nodes.append('')
            nodes = sorted(nodes, key=len, reverse=True)
            edges = [''.join(sorted(x['encoded_acids'])) for x in codons]
            g = nx.Graph()
            g.add_nodes_from(nodes)
            _add_graph_edges(g, nodes, edges, codons)
        else:
            g = self.graph

        for path in nx.all_shortest_paths(g, '', aa_str):
            codons_list = [
                g[v1][v2]['codon'] for v1, v2 in zip(path[:-1], path[1:])
            ]
            codon_comb = _combine_condons(codons_list)
            if codon_comb['score'] > best_codon_score:
                best_codon = codon_comb
                best_codon_score = codon_comb['score']

        return best_codon

    def _optimise_codons_greedy(self, amino_acids: tp.List[str],
                                codons: tp.List[tp.Dict]
                                ) -> tp.Optional[tp.Dict]:
        req_aas_set = set(amino_acids)

        for codons_list in _powerset(codons):
            codon_comb = _combine_condons(codons_list)
            if _is_codon_set_valid(req_aas_set, codon_comb['encoded_acids']):
                return codon_comb
        return None

    def _optimise_codons_exact(self, amino_acids: tp.List[str],
                               codons: tp.List[tp.Dict]
                               ) -> tp.Optional[tp.Dict]:
        best_codon = None
        best_codon_score = np.NINF
        best_codon_len = np.inf
        req_aas_set = set(amino_acids)

        for codons_list in _powerset(codons):
            codon_comb = _combine_condons(codons_list)
            if not _is_codon_set_valid(req_aas_set,
                                       codon_comb['encoded_acids']):
                continue
            if len(codon_comb['ambiguous_codons']) < best_codon_len:
                best_codon = codon_comb
                best_codon_len = len(codon_comb['ambiguous_codons'])
                best_codon_score = codon_comb['score']
                continue
            if codon_comb['score'] > best_codon_score:
                best_codon = codon_comb
                best_codon_len = len(codon_comb['ambiguous_codons'])
                best_codon_score = codon_comb['score']
                continue
        return best_codon

    def optimise_codons(self, amino_acids: tp.List[str],
                        mode='graph') -> tp.Optional[tp.Dict]:

        if mode == 'graph' and self.graph is not None:
            return self._optimise_codons_graph(amino_acids, None)

        if self.codon_dict is not None:
            codons = self._lookup_degenerate_codons(amino_acids)
        else:
            codons = self._gen_degenerate_codons(amino_acids)
            codons = _remove_duplicate_codons(codons)
        codons = sorted(codons,
                        key=lambda x: len(x['encoded_acids']),
                        reverse=True)

        if mode == 'graph':
            return self._optimise_codons_graph(amino_acids, codons)
        elif mode == 'greedy':
            return self._optimise_codons_greedy(amino_acids, codons)
        elif mode == 'exact':
            return self._optimise_codons_exact(amino_acids, codons)
        else:
            raise ValueError(
                'Unknown value %s for mode given.\
                    Allowed value are: graph, greedy, exact', mode)


P = tp.TypeVar('P')


def _powerset(p: tp.Sequence[P]) -> tp.Generator[tp.List[P], None, None]:
    for i in range(1, len(p) + 1):
        for item in combinations(p, i):
            yield list(item)


def _codon_to_codons_gen(codon_gen: tp.Generator[tp.Any, None, None]
                         ) -> tp.Generator[tp.Dict, None, None]:
    for codon in codon_gen:
        yield _codon_to_codons(codon)


def _codon_to_codons(codon: tp.Dict) -> tp.Dict:
    return {
        'ambiguous_codons': [codon['ambiguous_codon']],
        'amino_acids':
        codon['amino_acids'],
        'encoded_acids':
        set(amino_acid['amino_acid'] for amino_acid in codon['amino_acids']),
        'score':
        codon['score'],
    }


def _combine_condons(codons_list: tp.Sequence[tp.Dict]) -> tp.Dict:
    codons_dict = {
        'ambiguous_codons':
        list(chain(*(codons['ambiguous_codons'] for codons in codons_list))),
        'amino_acids':
        list(chain(*(codons['amino_acids'] for codons in codons_list))),
        'encoded_acids':
        set((aa for codons in codons_list for aa in codons['encoded_acids']))
    }
    codons_dict['score'] = _get_score(codons_dict['amino_acids'])
    return codons_dict


def _remove_duplicate_codons(codons_list: tp.Sequence[tp.Dict]
                             ) -> tp.List[tp.Dict]:
    codon_dict: tp.Dict[str, tp.Dict] = {}
    for codon in codons_list:
        try:
            if codon_dict[''.join(sorted(
                    codon['encoded_acids']))]['score'] < codon['score']:
                codon_dict[''.join(sorted(codon['encoded_acids']))] = codon
        except KeyError:
            codon_dict[''.join(sorted(codon['encoded_acids']))] = codon
    return list(codon_dict.values())


def _is_codon_set_valid(req_aas: set, codon_aas: set) -> bool:
    return len(req_aas.symmetric_difference(codon_aas)) == 0


def _add_graph_edges(g: nx.Graph, nodes: tp.List[str], edges: tp.List[str],
                     codons: tp.List[tp.Dict]):
    for v1, v2 in combinations(nodes, 2):
        diff = ''.join(sorted(set(v1) - set(v2)))
        try:
            idx = edges.index(diff)
            g.add_edge(v2, v1, codon=codons[idx])
        except ValueError:
            pass
