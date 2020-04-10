"""
CodonGenie (c) University of Manchester 2016

CodonGenie is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
"""
import itertools
from collections import defaultdict
from typing import Dict, Iterable, List, NamedTuple, Tuple

import Bio.Data.CodonTable as CodonTable
from synbiochem.utils import seq_utils
from synbiochem.utils.seq_utils import CodonOptimiser


class Codon(NamedTuple):
    codon: str
    probability: float
    cai: float


class AminoAcid(NamedTuple):
    type: int = -1
    amino_acid: str = "Stop"
    codons: Tuple[Codon, ...] = ()


class AmbigCodon(NamedTuple):
    ambiguous_codon: str
    ambiguous_codon_nucleotides: Tuple[str, ...]
    ambiguous_codon_expansion: Tuple[str, ...]
    amino_acids: Tuple[AminoAcid, ...]
    score: float = 0


class CodonSelector:
    """Class to optimise codon selection."""

    def __init__(self, table_id=1):
        self.__codon_to_aa = CodonTable.unambiguous_dna_by_id[table_id].forward_table
        self.__aa_to_codon = defaultdict(list)

        for codon, amino_acid in self.__codon_to_aa.items():
            self.__aa_to_codon[amino_acid].append(codon)

        self.__codon_opt = {}

    def optimise_codons(self, amino_acids, organism_id):
        """Optimises codon selection."""
        req_amino_acids = set(amino_acids.upper())

        codons = [seq_utils.CODONS[amino_acid] for amino_acid in req_amino_acids]

        results = [
            self.__analyse(combo, organism_id, req_amino_acids)
            for combo in itertools.product(*codons)
        ]

        return _format_results(results)

    def analyse_codon(self, ambig_codon, tax_id):
        """Analyses an ambiguous codon."""
        results = [[self.__analyse_ambig_codon(ambig_codon.upper(), tax_id)]]
        return _format_results(results)

    def __analyse(self, combo, tax_id, req_amino_acids):
        """Analyses a combo, returning nucleotides, ambiguous nucleotides,
        amino acids encodes, and number of variants."""
        transpose = [sorted(list(term)) for term in map(set, zip(*combo))]

        nucls = [["".join(sorted(list(set(pos))))] for pos in transpose[:2]] + [
            _optimise_pos_3(transpose[2])
        ]

        ambig_codons = [
            "".join([seq_utils.NUCL_CODES[term] for term in cdn])
            for cdn in itertools.product(*nucls)
        ]

        results = [
            self.__analyse_ambig_codon(ambig_codon, tax_id, req_amino_acids)
            for ambig_codon in ambig_codons
        ]

        return results

    def __analyse_ambig_codon(self, ambig_codon, tax_id, req_amino_acids=None):
        """Analyses a given ambiguous codon."""
        if req_amino_acids is None:
            req_amino_acids = []

        ambig_codon_nucls = [seq_utils.INV_NUCL_CODES[nucl] for nucl in ambig_codon]

        codons = ["".join(c) for c in itertools.product(*ambig_codon_nucls)]

        # amino_acids = defaultdict(dict)
        amino_acids = defaultdict(AminoAcid)

        for codon in codons:
            self.__analyse_codon(codon, tax_id, req_amino_acids, amino_acids)

        # amino_acids = [
        #     dict(val, **{"amino_acid": key})
        #     for key, val in sorted(
        #         amino_acids.items(), key=lambda x: (-x[1]["type"], x[0])
        #     )
        # ]

        # result = {
        #     "ambiguous_codon": ambig_codon,
        #     "ambiguous_codon_nucleotides": tuple(ambig_codon_nucls),
        #     "ambiguous_codon_expansion": tuple(codons),
        #     "amino_acids": amino_acids,
        # }

        result = AmbigCodon(
            ambiguous_codon=ambig_codon,
            ambiguous_codon_nucleotides=tuple(ambig_codon_nucls),
            ambiguous_codon_expansion=tuple(codons),
            amino_acids=tuple(amino_acids.values()),
        )

        if req_amino_acids:
            result = result._replace(score=_get_score(result.amino_acids))
            # result.update({"score": _get_score(amino_acids)})

        return result

    def __analyse_codon(
        self, codon, tax_id, req_amino_acids, amino_acids: Dict[str, AminoAcid]
    ):
        """Analyses a specific codon."""
        codon_opt = self.__get_codon_opt(tax_id)
        amino_acid = self.__codon_to_aa.get(codon, "Stop")
        typ = _get_amino_acid_type(amino_acid, req_amino_acids)
        aa = amino_acids[amino_acid]

        codon_tup = Codon(
            codon=codon,
            probability=codon_opt.get_codon_prob(codon),
            cai=codon_opt.get_cai(codon),
        )

        amino_acids[amino_acid] = aa._replace(
            type=typ, codons=aa.codons + (codon_tup,), amino_acid=amino_acid
        )

        # amino_acids[amino_acid]["type"] = typ

        # if "codons" not in amino_acids[amino_acid]:
        #     amino_acids[amino_acid]["codons"] = []
        #
        # amino_acids[amino_acid]["codons"].append(
        #     {
        #         "codon": codon,
        #         "probability": codon_opt.get_codon_prob(codon),
        #         "cai": codon_opt.get_cai(codon),
        #     }
        # )

    def __get_codon_opt(self, tax_id):
        """Gets the CodonOptimiser for the supplied taxonomy."""
        if tax_id not in self.__codon_opt:
            self.__codon_opt[tax_id] = CodonOptimiser(tax_id)

        return self.__codon_opt[tax_id]


def _optimise_pos_3(options):
    options = list({tuple(sorted(set(opt))) for opt in itertools.product(*options)})
    options.sort(key=len)
    return ["".join(opt) for opt in options]


def _get_amino_acid_type(amino_acid, req_amino_acids):
    """Gets amino acid type."""
    return -1 if amino_acid == "Stop" else (1 if amino_acid in req_amino_acids else 0)


def _get_score(amino_acids: Iterable[AminoAcid]):
    """Scores a given amino acids collection."""
    for vals in amino_acids:
        vals._replace(codons=tuple(sorted(vals.codons, key=lambda x: -x.cai)))
        # vals["codons"] = sorted(vals["codons"], key=lambda x: -x["cai"])

    # scores = [
    #     codon["cai"] if amino_acid["type"] == 1 else 0
    #     for amino_acid in amino_acids
    #     for codon in amino_acid["codons"]
    # ]
    scores = [
        codon.cai if amino_acid.type == 1 else 0
        for amino_acid in amino_acids
        for codon in amino_acid.codons
    ]

    return sum(scores) / float(len(scores))


def _format_results(results: List[List[AmbigCodon]]) -> List[AmbigCodon]:
    """Formats results."""
    return sorted(
        [codon for result in results for codon in result],
        key=lambda x: (len(x.ambiguous_codon_expansion), -x.score),
    )
    # return sorted(
    #     [codon for result in results for codon in result],
    #     key=lambda x: (
    #         len(x["ambiguous_codon_expansion"]),
    #         -x["score"] if "score" in x else 0,
    #     ),
    # )
