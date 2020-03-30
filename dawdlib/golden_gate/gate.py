from typing import NamedTuple, Optional, Tuple, Union

from pandas import DataFrame


class SynMut(NamedTuple):
    idx: int
    aa: str
    codon: str


class Gate(NamedTuple):
    idx: Optional[int]
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
        return (
            self.idx == other.idx
            and self.bps == other.bps
            and self.req_primer == other.req_primer
            and self.syn_mut == other.syn_mut
        )

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

    def span(self) -> Union[Tuple[int, int], Tuple[None, None]]:
        try:
            return self.idx, self.idx + 3
        except TypeError:
            return None, None


class GateSet(NamedTuple):
    no_oligos: str
    idx: int
    table: DataFrame
