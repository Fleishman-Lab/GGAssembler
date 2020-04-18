from heapq import heappop, heappush
from typing import Dict, Iterable, List, Set, Tuple

from dawdlib.degenerate_dna.structs import PosCodon

AADict = Dict[str, List[Tuple[int, PosCodon]]]
PosDict = Dict[PosCodon, List[str]]
Row = PosCodon
Col = str
Cols = Dict[str, List[PosCodon]]


class AlgXSolver:
    max_size: int

    @classmethod
    def solve(
        cls, amino_acids: Set[str], codons: List[PosCodon]
    ) -> Iterable[List[PosCodon]]:
        """
        Solve the exact cover problem, using Algorith X.
        Args:
            amino_acids:
            codons:

        Returns:
        """
        aa_d, cdn_d = cls._prep_data(amino_acids, codons)
        cls.max_size = len(amino_acids)
        return cls._solve(aa_d, cdn_d)

    @staticmethod
    def _prep_data(
        amino_acids: Set[str], codons: List[PosCodon]
    ) -> Tuple[AADict, PosDict]:
        amb_cdn: Tuple[str, ...]
        pos_cdn: PosCodon
        cdn_d: PosDict = dict((p, list(p.encoded_acids)) for p in codons)
        aa_d: AADict = dict((aa, []) for aa in amino_acids)
        for pos_cdn, enc_aas in cdn_d.items():
            for aa in enc_aas:
                heappush(aa_d[aa], (-len(pos_cdn.encoded_acids), pos_cdn))
                # aa_d[aa].add(pos_cdn)
        return aa_d, cdn_d

    @classmethod
    def _solve(
        cls, aa_d: AADict, cdn_d: PosDict, solution: List[PosCodon] = None,
    ) -> Iterable[List[PosCodon]]:
        c: Col
        r: Row
        if solution is None:
            solution = []
        if len(solution) > cls.max_size:
            return
        elif not aa_d:
            cls.max_size = len(solution)
            yield solution
        else:
            c = min(aa_d, key=lambda k: len(aa_d[k]))
            for _, r in aa_d[c]:
                solution.append(r)
                cols = cls._select(aa_d, cdn_d, r)
                for s in cls._solve(aa_d, cdn_d, solution):
                    cls.max_size = len(s)
                    yield s
                cls._deselect(aa_d, cdn_d, r, cols)
                solution.pop()

    @staticmethod
    def _select(aa_d: AADict, cdn_d: PosDict, r: Row) -> Cols:
        cols: Cols = {}
        for j in cdn_d[r]:
            for idx, i in aa_d[j]:
                for k in cdn_d[i]:
                    if k != j:
                        aa_d[k].remove((idx, i))
            cols[j] = aa_d.pop(j)
        return cols

    @staticmethod
    def _deselect(aa_d: AADict, cdn_d: PosDict, r: Row, cols: Cols):
        for j in reversed(cdn_d[r]):
            aa_d[j] = cols.pop(j)
            for idx, i in aa_d[j]:
                for k in cdn_d[i]:
                    if k != j:
                        heappush(aa_d[k], (-len(i.encoded_acids), i))
                        # aa_d[k].he(i)
