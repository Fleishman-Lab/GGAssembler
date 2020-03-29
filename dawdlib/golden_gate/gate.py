from typing import NamedTuple, Tuple

from pandas import DataFrame


class SynMut(NamedTuple):
    idx: int
    aa: str
    codon: str


class Gate(NamedTuple):
    idx: int
    bps: str
    req_primer: bool = False
    syn_mut: Tuple[SynMut, ...] = ()

    def __lt__(self, other: "Gate") -> bool:
        if not isinstance(other, Gate):
            return False
        return self.idx < other.idx

    def __le__(self, other: "Gate") -> bool:
        if not isinstance(other, Gate):
            return False
        return self.idx <= other.idx

    def __eq__(self, other: "Gate") -> bool:
        if not isinstance(other, Gate):
            return False
        return self.idx == other.idx

    def __ne__(self, other: "Gate") -> bool:
        if not isinstance(other, Gate):
            return False

        return self.idx != other.idx

    def __gt__(self, other: "Gate") -> bool:
        if not isinstance(other, Gate):
            return False

        return self.idx > other.idx

    def __ge__(self, other: "Gate") -> bool:
        if not isinstance(other, Gate):
            return False

        return self.idx >= other.idx

    def __sub__(self, other: "Gate") -> 0:
        if not isinstance(other, Gate):
            return False
        if self.idx == other.idx:
            return 0
        if self.idx < other.idx:
            return self.idx - other.idx - 4
        if self.idx > other.idx:
            return self.idx + 4 - other.idx

    def overlap(self, other: "Gate") -> bool:
        if not isinstance(other, Gate):
            return False
        if self.idx == other.idx:
            return True
        if self.idx < other.idx:
            return other.idx <= self.idx + 3
        return self.idx <= other.idx + 3

    def span(self) -> Tuple[int, int]:
        return self.idx, self.idx + 4


class GateSet(NamedTuple):
    no_oligos: str
    idx: int
    table: DataFrame
