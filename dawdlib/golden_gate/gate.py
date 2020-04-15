from typing import NamedTuple, Tuple

from pandas import DataFrame


class SynMut(NamedTuple):
    idx: int
    aa: str
    codon: str


class Gate(NamedTuple):
    idx: int
    bps: str
    src_or_target: bool = False
    req_primer: bool = False
    syn_mut: Tuple[SynMut, ...] = ()

    def __sub__(self, other: "Gate") -> int:
        # 4 is added to account for gate length (4)
        # plus the fact we need to include the entire gate
        if self.idx == other.idx:
            return 0
        if self < other:
            return self.idx - other.idx - 4
        # only option left is self.idx > other.idx:
        return self.idx + 4 - other.idx

    def overlap(self, other: "Gate") -> bool:
        if self == other:
            return True
        if self < other:
            return other.idx <= self.idx + 3
        return self.idx <= other.idx + 3

    def span(self) -> Tuple[int, int]:
        if self.src_or_target:
            return self.idx, self.idx
        return self.idx, self.idx + 3


class GateSet(NamedTuple):
    no_oligos: str
    idx: int
    table: DataFrame
