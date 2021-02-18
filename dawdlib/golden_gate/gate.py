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
    gatelength: int = 4

    def __sub__(self, other: "Gate") -> int:
        # 4 is added to account for gate length (4)
        # plus the fact we need to include the entire gate
        diff = self.idx - other.idx
        if diff:
            return abs(diff) + gatelength
        return diff

    def overlap(self, other: "Gate") -> bool:
        if isinstance(other, PseudoGate):
            return False
        if self == other:
            return True
        if self < other:
            return other.idx <= self.idx + gatelength - 1
        return self.idx <= other.idx + gatelength - 1

    def span(self) -> Tuple[int, int]:
        if self.src_or_target:
            return self.idx, self.idx
        return self.idx, self.idx + gatelength - 1


class PseudoGate(Gate):
    def __sub__(self, other: "PseudoGate") -> int:
        diff = super().__sub__(other)
        if diff:
            return diff - 1
        return diff

    def overlap(self, other: "PseudoGate") -> bool:
        return False

    def span(self) -> Tuple[int, int]:
        return self.idx, self.idx


class GateSet(NamedTuple):
    no_oligos: str
    idx: int
    table: DataFrame
