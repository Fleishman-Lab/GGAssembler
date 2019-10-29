import pickle
import typing as tp

from dawdlib.degenerate_dna.codon_selector import _generate_codons_graph


def write_codons_graph(
    amino_acids: tp.List[str], codons: tp.List[tp.Dict], filename: str
) -> None:
    g = _generate_codons_graph(amino_acids, codons)
    with open(filename, "wb") as outfile:
        pickle.dump(g, outfile)
