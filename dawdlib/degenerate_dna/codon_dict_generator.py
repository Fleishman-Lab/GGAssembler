import pickle
import typing as tp

from dawdlib.degenerate_dna.codon_selector import (
    CodonSelector,
    _remove_duplicate_codons,
)


def generate_codon_dict(organism_id: str, amino_acids: tp.List[str]) -> tp.Dict:
    cdn_sel = CodonSelector(organism_id)
    codons = cdn_sel._gen_degenerate_codons(amino_acids)
    codons = _remove_duplicate_codons(codons)
    return dict(("".join(sorted(codon.encoded_acids)), codon) for codon in codons)


def write_codon_dict(
    organism_id: str, amino_acids: tp.List[str], filename: str
) -> None:
    codon_dict = generate_codon_dict(organism_id, amino_acids)
    with open(filename, "wb") as outf:
        pickle.dump(codon_dict, outf)
