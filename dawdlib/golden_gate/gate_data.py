import os
from itertools import combinations
from typing import List

import pandas as pd


class GGData:

    default_df = "%s/resources/FileS01_T4_18h_37C.csv" % os.path.dirname(__file__)

    def __init__(self, init_df=True) -> None:
        self.lig_df: pd.DataFrame = None
        if init_df:
            self.lig_df = self._parse_ligation_df()

    def _load_ligation_df(self, df_path: str, *args, **kwargs) -> pd.DataFrame:
        return pd.read_csv(df_path, *args, **kwargs)

    def _parse_ligation_df(self):
        return self._load_ligation_df(self.default_df, index_col=0)

    def gates_scores(self, gate1: str, gate2: str) -> int:
        score1 = self.lig_df.loc[gate1, gate2]
        score2 = self.lig_df.loc[gate2, gate1]
        return int((score1 + score2) / 2)

    def gates_all_scores(self, gate1: str, gate2: str) -> int:
        score = 0
        score += self.gates_scores(gate1, gate2)
        score += self.gates_scores(gate1, reverse_complement(gate2))
        score += self.gates_scores(gate2, reverse_complement(gate1))
        score += self.gates_scores(reverse_complement(gate1), reverse_complement(gate2))
        return score

    def score_gate_clique(self, gates_clq: List[str]) -> int:
        score = 0
        for gate1 in gates_clq:
            for gate2 in gates_clq:
                if gate1 != gate2:
                    score += self.gates_all_scores(gate1, gate2)
        return score

    def score_on_gate_clique(self, gates: List[str]) -> int:
        return sum([self.gates_scores(g, reverse_complement(g)) for g in gates])

    def is_one_over(self, gates: List[str], threshold: int = 1000) -> bool:
        """used for screening lists of gates where at least one off target gate pair scores over threshold

        Parameters
        ----------
        gates : List[str]
            a list of gates (clique) to be testes
        threshold : int, optional
            threshold to consider a gate pair as off-target binders, by default 1000

        Returns
        -------
        bool
            if True, the gates list has an off-target pair that scores over threshold
        """
        for gate1 in gates:
            for gate2 in gates:
                if gate1 != gate2:
                    if self.gates_scores(gate1, gate2) > threshold:
                        return True
        return False

    def is_one_on_target_under(self, gates: List[str], threshold: int = 2000) -> bool:
        """used for screening lists of gates where at least one gate binds itself under threshold

        Parameters
        ----------
        gates : List[str]
            list of gates, clique
        threshold : int, optional
            threshold to consider as binders, by default 2000

        Returns
        -------
        bool
            whether the gates list has at least one gate that does not bind itself well
        """
        for gate in gates:
            if self.gates_scores(gate, reverse_complement(gate)) < threshold:
                return True
        return False

    def score_gate1_rc_gate2(self, gate1: str, gate2: str) -> int:
        return self.gates_scores(gate1, reverse_complement(gate2))

    def self_binding_scores(self) -> List[int]:
        return [self.score_gate1_rc_gate2(g, g) for g in self.lig_df.columns]

    def filter_self_binding_gates(self, threshold: int = 2000) -> List[str]:
        return [
            g
            for g, score in zip(self.lig_df.columns, self.self_binding_scores())
            if score > threshold and g != reverse_complement(g)
        ]

    def gate_set_has_off_target(self, gates: List[str], threshold: int = 1000) -> bool:
        for gate1, gate2 in combinations(gates, 2):
            if self.gates_all_scores(gate1, gate2) > threshold:
                return True
        return False


def reverse_complement(seq) -> str:
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join([complement[a] for a in seq[::-1]])
