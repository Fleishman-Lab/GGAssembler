import typing as tp
from enum import Enum

DEG_NUCL_CODES = [
    "A",
    "C",
    "G",
    "T",
    "R",
    "Y",
    "S",
    "W",
    "K",
    "M",
    "B",
    "D",
    "H",
    "V",
    "N",
]


class Codon(tp.NamedTuple):
    codon: str
    probability: float
    cai: float


class AminoAcid(tp.NamedTuple):
    type: int = -1
    amino_acid: str = "Stop"
    codons: tp.Tuple[Codon, ...] = ()


class AmbigCodon(tp.NamedTuple):
    ambiguous_codon: str
    ambiguous_codon_nucleotides: tp.Tuple[str, ...]
    ambiguous_codon_expansion: tp.Tuple[str, ...]
    amino_acids: tp.Tuple[AminoAcid, ...]
    score: float = 0


class PosCodon(tp.NamedTuple):
    ambiguous_codons: tp.Tuple[str, ...] = ()
    amino_acids: tp.Tuple[AminoAcid, ...] = ()
    encoded_acids: tp.FrozenSet[str] = frozenset()
    score: float = -1


class TableColNames(Enum):
    AA_POS = "AA_POS"
    DNA_POS = "DNA_POS"
    AMBIGUOUS_CODONS = "AMBIGUOUS_CODONS"
    ENCODED_AAS = "ENCODED_AAS"
    ENCODED_COUNT = "ENCODED_COUNT"
