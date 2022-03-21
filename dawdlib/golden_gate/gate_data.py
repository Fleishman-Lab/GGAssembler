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
from typing import DefaultDict, Dict, FrozenSet, Iterable, List, NamedTuple, Tuple

import numpy as np
import pandas as pd

import dawdlib.golden_gate.resources as gg_resources


class FidelityResults(NamedTuple):
    NEB: float
    FWD: float
    REV: float


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
        temperature: int = 25,
        hours: int = 18,
        min_efficiency: float = MIN_EFFICIENCY,
        min_fidelity: float = MIN_FIDELITY,
    ) -> None:
        self.lig_df: pd.DataFrame

        self._min_efficiency: float = min_efficiency
        self._min_fidelity: float = min_fidelity
        self._efficient_overhangs: List[str] = []
        self.lig_df: pd.DataFrame = pd.DataFrame()
        self.fwd_efficiency: pd.DataFrame = pd.DataFrame()
        self.rev_efficiency: pd.DataFrame = pd.DataFrame()
        self.fwd_fidelity: pd.DataFrame = pd.DataFrame()
        self.rev_fidelity: pd.DataFrame = pd.DataFrame()
        self.initialized = False
        try:
            self.default_df = self.df_files[(hours, temperature)]
            self.init()
        except KeyError:
            raise ValueError(
                "No data was found for the combination of temperature %d C and time %d H."
                % (temperature, hours)
            )

        self.score_dict: Dict[FrozenSet, int] = {}
        self.all_score_dict: Dict[FrozenSet, int] = {}
        self._rev_dict: Dict[str, str] = {}

    def init(self) -> None:
        self.lig_df = self._parse_ligation_df().sort_index(axis=0)
        fwd_efficiency = self.lig_df.sum(axis=1).sort_index(axis=0)
        rev_efficiency = self.lig_df.sum(axis=0).sort_index(axis=0)
        self.fwd_efficiency = fwd_efficiency / fwd_efficiency.max()
        self.rev_efficiency = rev_efficiency / rev_efficiency.max()
        self.fwd_fidelity = self.lig_df
        self.rev_fidelity = self.lig_df
        self.initialized = True

    def set_default_df(self, csv: str) -> None:
        if os.path.exists(csv):
            self.default_df = csv
        else:
            self.default_df = gg_resources.ligation_data[csv]

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

    def _parse_ligation_df(self) -> pd.DataFrame:
        df = self._load_ligation_df(self.default_df, index_col=0)
        # According to NEB we need to transform values so the sum os 100K
        return (df * 1e5) / df.values.sum()

    def filter_self_binding_gates(self, filter_gc: bool = True) -> List[str]:
        if not self._efficient_overhangs:
            overhangs = self.fwd_efficiency
            if filter_gc:
                overhangs = self.fwd_efficiency.filter(regex=r"[A|T]")
            revs = sorted(overhangs.index.map(reverse_complement))
            overhangs = overhangs[
                (overhangs > self.min_efficiency).values
                & (self.rev_efficiency.loc[revs] > self.min_efficiency).values
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
        self.fwd_fidelity = sub_lig_df.sort_index(axis=0)
        self.rev_fidelity = sub_lig_df.T.sort_index(axis=0)

    def restriction_edges(self, fwd_overhangs: List[str]) -> Iterable[Tuple[str, str]]:
        allowed_mismatch = 10 * (1 - self.min_fidelity)
        rev_overhangs = list(map(reverse_complement, fwd_overhangs))
        overhangs = fwd_overhangs + rev_overhangs
        for over in overhangs:
            for fid_df in [self.fwd_fidelity, self.rev_fidelity]:
                fidelity = fid_df[over]
                for ligating_over in fidelity[fidelity > allowed_mismatch].index:
                    yield over, ligating_over

    def overhangs_fidelity(self, *args) -> Iterable[float]:
        fwds = list(map(str, args))
        revs = list(map(reverse_complement, fwds))
        overhangs = fwds + revs
        for over, rev in zip(fwds, revs):
            if fwds.count(over) > 1:
                yield 0.0
            elif rev in fwds:
                yield 0.0
            else:
                yield min(
                    self.fwd_fidelity.loc[over, rev]
                    / self.fwd_fidelity.loc[over, overhangs].sum(),
                    self.fwd_fidelity.loc[rev, over]
                    / self.fwd_fidelity.loc[rev, overhangs].sum(),
                    self.rev_fidelity.loc[over, rev]
                    / self.rev_fidelity.loc[over, overhangs].sum(),
                    self.rev_fidelity.loc[rev, over]
                    / self.rev_fidelity.loc[rev, overhangs].sum(),
                )

    def reaction_fidelity(self, *args) -> Tuple[float, float, float]:
        fwds = sorted(map(str, args))
        revs = list(map(reverse_complement, fwds))
        fwds_set = set(fwds)
        revs_set = set(revs)
        if len(fwds) != len(fwds_set):
            raise ValueError("At least one overhangs appears more than once.")
        if fwds_set.intersection(revs_set):
            raise ValueError(
                "At least one overhang and it's reverse complement are present."
            )

        fwd_fidelity = self._directional_fidelity(fwds, revs, True)
        rev_fidelity = self._directional_fidelity(revs, fwds, False)
        reaction_fidelity = self._bidirectional_fidelity(fwds, revs)

        return FidelityResults(fwd_fidelity, rev_fidelity, reaction_fidelity)

    def _directional_fidelity(self, fwds: List[str], revs: List[str], fwd=True):
        fidelity_df = self.fwd_fidelity if fwd else self.rev_fidelity
        overhangs = fwds + revs
        all_fid = (fidelity_df.loc[fwds, overhangs]).sum(axis=1, skipna=True)
        complement = fidelity_df.loc[fwds, revs].values.diagonal()
        return complement.divide(all_fid).prod()

    def _bidirectional_fidelity(self, fwds: List[str], revs: List[str]):
        fidelity_df = self.fwd_fidelity
        over_pairs = list(zip(fwds, revs))
        overhangs = fwds + revs
        pair_ligations = np.array(
            [
                (fidelity_df.loc[over_pair, overhangs]).sum(skipna=True)
                for over_pair in over_pairs
            ]
        )
        complement = fidelity_df.loc[fwds, revs].values.diagonal()
        return complement.divide(pair_ligations - complement).prod()

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
                    cuts_min_fidelity, self.reaction_fidelity(overhangs[a_cut:b_cut])[0]
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
