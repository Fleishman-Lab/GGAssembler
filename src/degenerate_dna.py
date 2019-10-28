import typing as tp
from functools import reduce
from itertools import chain, combinations

import networkx as nx
import numpy as np
from codon_utils import CodonSelector as UtilsCodonSelector
from codon_utils import _get_score


class CodonSelector(object):
    def __init__(self, organism_id: str):
        self.organism_id = organism_id
        self.cod_sel = UtilsCodonSelector()

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

    def optimise_codons(self, amino_acids: tp.List[str],
                        exact=False) -> tp.Optional[tp.Dict]:
        '''Build a graph of possible degenerate codons
        that only encode the required given codons.
        Find shortest weighted path to get best codon combinations.

        Error is defind as:
            the number of aas provided by the codon
            minus the number of required codons to the power of 3.

            e = (|amino_acids| - |codon_aas|)^3

        An edge is connected only if the the codons create the entire
        set of required AAs.
        '''

        best_codon = None
        best_codon_score = np.NINF
        # req_aas_set = set(amino_acids)
        aa_str = ''.join(sorted(amino_acids))

        # codons = list(
        #     chain(*(_codon_to_codons_gen(
        #         self.select_codons(''.join(amino_acid_set)))
        #             for amino_acid_set in _powerset(amino_acids))))
        codons = self._gen_degenerate_codons(amino_acids)
        codons = _remove_duplicate_codons(codons)
        codons = sorted(codons,
                        key=lambda x: len(x['encoded_acids']),
                        reverse=True)

        nodes = [''.join(sorted(x)) for x in _powerset(amino_acids)]
        nodes.append('')
        nodes = sorted(nodes, key=len, reverse=True)
        edges = [''.join(sorted(x['encoded_acids'])) for x in codons]
        g = nx.Graph()
        g.add_nodes_from(nodes)
        _add_graph_edges(g, nodes, edges, codons)

        for path in nx.all_shortest_paths(g, '', aa_str):
            codons_list = [g[v1][v2]['codon'] for v1, v2 in zip(path[:-1], path[1:])]
            codon_comb = _combine_condons(codons_list)
            if codon_comb['score'] > best_codon_score:
                best_codon = codon_comb
                best_codon_score = codon_comb['score']

        # if exact:
        #     for codons_list in _powerset(codons):
        #         codon_comb = _combine_condons(codons_list)
        #         if not _is_codon_set_valid(req_aas_set,
        #                                    codon_comb['encoded_acids']):
        #             continue
        #         if len(codon_comb['ambiguous_codons']) < best_codon_len:
        #             best_codon = codon_comb
        #             best_codon_len = len(codon_comb['ambiguous_codons'])
        #             best_codon_score = codon_comb['score']
        #             continue
        #         if codon_comb['score'] > best_codon_score:
        #             best_codon = codon_comb
        #             best_codon_len = len(codon_comb['ambiguous_codons'])
        #             best_codon_score = codon_comb['score']
        #             continue
        # else:
        #     for codons_list in _powerset(codons):
        #         codon_comb = _combine_condons(codons_list)
        #         if _is_codon_set_valid(req_aas_set,
        #                                codon_comb['encoded_acids']) and codon_comb['score'] == 1:
        #             return codon_comb

        return best_codon


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
    codon_dict = {}
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
