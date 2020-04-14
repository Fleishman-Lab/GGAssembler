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

    # Since gate is based on a tuple and only a single gate be at the a specific
    # index tuple comparison is enough to check all these cases.
    # def __lt__(self, other: Tuple) -> Union[bool, "NotImplemented"]:
    #     if not isinstance(other, Gate):
    #         return NotImplemented
    #     return self.idx < other.idx
    #
    # def __le__(self, other: Tuple) -> Union[bool, "NotImplemented"]:
    #     if not isinstance(other, Gate):
    #         return NotImplemented
    #     return self.idx <= other.idx
    #
    # def __eq__(self, other: object) -> Union[bool, "NotImplemented"]:
    #     if not isinstance(other, Gate):
    #         return NotImplemented
    #     return (
    #         self.idx == other.idx
    #         and self.bps == other.bps
    #         and self.req_primer == other.req_primer
    #         and self.syn_mut == other.syn_mut
    #     )
    #
    # def __ne__(self, other: object) -> Union[bool, "NotImplemented"]:
    #     if not isinstance(other, Gate):
    #         return NotImplemented
    #
    #     return self.idx != other.idx
    #
    # def __gt__(self, other: Tuple) -> Union[bool, "NotImplemented"]:
    #     if not isinstance(other, Gate):
    #         return NotImplemented
    #
    #     return self.idx > other.idx
    #
    # def __ge__(self, other: Tuple) -> Union[bool, "NotImplemented"]:
    #     if not isinstance(other, Gate):
    #         return NotImplemented
    #
    #     return self.idx >= other.idx

    def __sub__(self, other: "Gate") -> int:
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
