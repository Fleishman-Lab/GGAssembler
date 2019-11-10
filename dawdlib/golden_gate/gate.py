from typing import NamedTuple
from pandas import DataFrame


class Gate(NamedTuple):
    idx: int
    bps: str


class GateSet(NamedTuple):
    no_oligos: str
    idx: int
    table: DataFrame
