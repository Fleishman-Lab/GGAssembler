"""
References:
Potapov, V., et al. (2018).
    Comprehensive Profiling of Four Base Overhang Ligation Fidelity by T4 DNA Ligase and Application to DNA Assembly.
    ACS Synth. Biol. 7, 2665–2674.
"""

import math
import os
from collections import defaultdict
from functools import lru_cache
from typing import DefaultDict, Dict, FrozenSet, Iterable, List, Tuple

import pandas as pd


class GGData:

    df_files = {
        (18, 37): "%s/resources/FileS_T4_18h_37C.csv" % os.path.dirname(__file__),
        (18, 25): "%s/resources/FileS_T4_18h_25C.csv" % os.path.dirname(__file__),
        (1, 25): "%s/resources/FileS_T4_01h_25C.csv" % os.path.dirname(__file__),
    }

    MIN_EFFICIENCY = 0.25
    MIN_FIDELITY = 0.9

    def __init__(
        self,
        init_df=True,
        temperature: int = 37,
        hours: int = 18,
        min_efficiency: float = MIN_EFFICIENCY,
        min_fidelity: float = MIN_FIDELITY,
    ) -> None:
        self.lig_df: pd.DataFrame
        # Efficiency is proportional to the "best" overhang in the data.
        # It is the proportion of an overhang from the "best" overhang
        self.efficiency: pd.Series
        # Fidelity is the probability of an overhang to connect to it's reverse complement
        self.fidelity: pd.DataFrame

        self._min_efficiency: float = min_efficiency
        self._min_fidelity: float = min_fidelity
        self._efficient_overhangs: List[str] = []

        try:
            self.default_df = self.df_files[(hours, temperature)]
        except KeyError:
            raise ValueError(
                "No data was found for the combination of temperature %d C and time %d H."
                % (temperature, hours)
            )

        self.score_dict: Dict[FrozenSet, int] = {}
        self.all_score_dict: Dict[FrozenSet, int] = {}
        self._rev_dict: Dict[str, str] = {}

        if init_df:
            self.lig_df = self._parse_ligation_df()
            efficiency = self.lig_df.sum(axis=0)
            self.efficiency = efficiency / efficiency.max()
            self.fidelity = self.lig_df.divide(efficiency)

    def get_efficiency(self) -> float:
        return self._min_efficiency

    def set_efficiency(self, min_efficiency: float):
        self._min_efficiency = min_efficiency

    min_efficiency = property(get_efficiency, set_efficiency)

    def get_fidelity(self) -> float:
        return self._min_fidelity

    def set_fidelity(self, min_fidelity: float):
        self._min_fidelity = min_fidelity

    min_fidelity = property(get_fidelity, set_fidelity)

    @staticmethod
    def _load_ligation_df(df_path: str, *args, **kwargs) -> pd.DataFrame:
        return pd.read_csv(df_path, *args, **kwargs)

    def _parse_ligation_df(self):
        return self._load_ligation_df(self.default_df, index_col=0)

    def filter_self_binding_gates(self, filter_gc: bool = True) -> List[str]:
        if not self._efficient_overhangs:
            overhangs = self.efficiency
            if filter_gc:
                overhangs = self.efficiency.filter(regex=r"[A|T]")
            revs = overhangs.index.map(reverse_complement)
            overhangs = overhangs[
                (overhangs > self.min_efficiency).values
                & (overhangs.loc[revs] > self.min_efficiency).values
            ]
            self._efficient_overhangs = list(overhangs.index)
            self._update_fidelity()
        return self._efficient_overhangs

    def _update_fidelity(self):
        assert (
            self._efficient_overhangs
        ), "filter_self_binding_gates has not been called yet, cannot run!"
        sub_lig_df = self.lig_df.loc[
            self._efficient_overhangs, self._efficient_overhangs
        ]
        efficiency = sub_lig_df.sum(axis=0)
        self.fidelity = sub_lig_df.divide(efficiency)

    def restriction_edges(self, overhangs: List[str]) -> Iterable[Tuple[str, str]]:
        for over in overhangs:
            rev = reverse_complement(over)
            node_ser = self.fidelity.loc[:, over]
            yield over, rev
            if node_ser.loc[rev] >= self.min_fidelity or math.isclose(
                node_ser.loc[rev], self.min_fidelity, abs_tol=1e-2
            ):
                continue
            tot = 0
            cutoff = 1 - node_ser.loc[rev] / self.min_fidelity
            srtd_node_ser = node_ser.loc[node_ser.index != rev].sort_values(
                ascending=False
            )
            for indx, val in srtd_node_ser.iteritems():
                if val < 1e-2:
                    break
                tot += val
                yield over, indx
                if tot >= cutoff or math.isclose(cutoff, tot, abs_tol=1e-2):
                    break

    @lru_cache(maxsize=None)
    def overhangs_fidelity(self, *args) -> Iterable[float]:
        overhangs = list(map(str, args))
        revs = list(map(reverse_complement, overhangs))
        for over, rev in zip(overhangs, revs):
            if overhangs.count(over) > 1:
                yield 0.0
            elif rev in overhangs:
                yield 0.0
            else:
                yield min(
                    self.fidelity.loc[rev, over]
                    / self.fidelity.loc[revs + overhangs, over].sum(),
                    self.fidelity.loc[over, rev]
                    / self.fidelity.loc[revs + overhangs, rev].sum(),
                )

    def create_subreactions(
        self, overhangs: List[str], max_reactions: int
    ) -> DefaultDict[int, Tuple[float, FrozenSet[int]]]:
        best_reactions: DefaultDict[int, Tuple[float, FrozenSet[int]]] = defaultdict(
            lambda: (0, frozenset())
        )
        for cuts in nmulticombinations(len(overhangs), max_reactions):
            cuts_min_fidelity = 1
            for a_cut, b_cut in zip(cuts[:-1], cuts[1:]):
                cuts_min_fidelity = min(
                    cuts_min_fidelity,
                    min(self.overhangs_fidelity(overhangs[a_cut:b_cut])),
                )
                if cuts_min_fidelity < best_reactions[len(cuts)][0]:
                    continue
            best_reactions[len(cuts)] = (cuts_min_fidelity, cuts)
        return best_reactions


@lru_cache(maxsize=256)
def reverse_complement(seq) -> str:
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join([complement[a] for a in seq[::-1]])


def nmulticombinations(
    balls: int, bins: int, partition: Tuple[int, ...] = (0,)
) -> Iterable[FrozenSet[int]]:
    if bins <= 0:
        raise ValueError()
    if bins == 1:
        yield sorted(frozenset(partition + (balls,)))
    else:
        for b in range(partition[-1], balls):
            yield from nmulticombinations(balls, bins - 1, partition + (b,))
