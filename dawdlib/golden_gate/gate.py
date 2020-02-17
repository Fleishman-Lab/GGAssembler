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


class GateSet(NamedTuple):
    no_oligos: str
    idx: int
    table: DataFrame
