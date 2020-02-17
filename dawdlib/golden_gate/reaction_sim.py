from typing import Generator, Iterator, NamedTuple, Tuple

from Bio.Restriction.Restriction import Ov5
from Bio.Seq import Seq


class OverHang(NamedTuple):
    fprime: int = 0
    tprime: int = 0


class OverHangDNA(NamedTuple):
    fprime: str = ""
    tprime: str = ""


class SDNASection(NamedTuple):
    dna: str = ""
    start: int = 0
    end: int = 0
    hang: OverHang = OverHang()
    hang_dna: OverHangDNA = OverHangDNA()


class DDNASection(NamedTuple):
    fprime: SDNASection = SDNASection()
    tprime: SDNASection = SDNASection()


def get_ddna_sections(
    enzyme: Ov5, dna_segs: Iterator[Tuple[str, str]], start: int
) -> Generator[DDNASection, None, None]:
    """

    Args:
        enzyme (Ov5): Restriction enzyme which was used to digest the DNA.
        dna_segs (Iterator[Tuple[str, str]]): An iterator over double stranded DNA segments
                                              (after digestion by restrction enzyme).
        start (int): An index pointing to the start of the DNA segment or the location of the gate.

    Yields:
        DDNASection: The double stranded DNA segment.
    """
    prev_ddna = DDNASection()
    prev_ddna.five_three.end = start
    for fe, te in dna_segs:
        f_overhang, t_overhang = find_overhang(
            len(fe), len(te), prev_ddna.fprime.hang.tprime, prev_ddna.tprime.hang.fprime
        )
        f_dna_overhang, t_dna_overhang = find_overhand_dna(
            fe, te, f_overhang, t_overhang
        )
        prev_ddna = DDNASection(
            SDNASection(
                fe,
                prev_ddna.fprime.end,
                prev_ddna.fprime + len(fe),
                f_overhang,
                f_dna_overhang,
            ),
            SDNASection(dna=te, hang=t_overhang, hang_dna=t_dna_overhang),
        )
        # Skip yielding sections which have the enzyme recognition site
        if enzyme.search(fe) or enzyme.search(te):
            continue
        yield prev_ddna


def enzyme_dna_segments(enzyme: Ov5, dna: str) -> Iterator[Tuple[str, str]]:
    """

    Args:
        enzyme (Ov5): Restriction enzyme used to cut the dna
        dna (str): 5' DNA string

    Returns:
        Iterator[Tuple[str, str]]: An iterator over the double stranded DNA segments
    """
    seq = Seq(dna)
    f2t = map(str, enzyme.catalyse(seq))
    t2f = reversed(
        list(
            map(
                lambda x: "".join(reversed(x)),
                enzyme.catalyse(seq.reverse_complement()),
            )
        )
    )
    return iter(zip(f2t, t2f))


def find_overhang(
    flen: int, tlen: int, pft: int, ptf: int
) -> Tuple[OverHang, OverHang]:
    """
    Find DNA hanging ends.

    Takes care of six different cases:
        1. ----------
           ---------------
           flen: 0, tlen: 0 -> flen: 0, tlen: 5
        2. ---------------
                ----------
           flen: 0, tlen: 5 -> flen: 0, tlen: 0
        3. ---------------
                ---------------
           flen: 0, tlen: 5 -> flen: 0, tlen: 5
        And the reverse of these cases
    Args:
        flen (int): Length of DNA string of 5' to 3'
        tlen (int): Length of DNA string of 3' to 5'
        pft (int): previous 5' segment 3' side hanging length
        ptf (int): previous 3' segment 5' side hanging length

    Returns:
        Tuple[OverHang, OverHang]: A tuple with the lengths of the hanging DNA of both stands (5' and 3').
    """

    cft = flen + pft - tlen
    cft = cft if cft >= 0 else 0
    ctf = tlen + ptf - flen
    ctf = ctf if ctf >= 0 else 0
    return OverHang(ptf, cft), OverHang(pft, ctf)


def find_overhand_dna(
    fdna: str, tdna: str, f_overhang: OverHang, t_overhang: OverHang
) -> Tuple[OverHangDNA, OverHangDNA]:
    """
    See ~find_overhang for cases and logic.

    Args:
        fdna (str): The 5' DNA string
        tdna (str): The 3' DNA string:
        f_overhang (OverHang): The 5' strand overhangs
        t_overhang (OverHang): The 5' strand overhangs

    Returns:
        Tuple[OverHangDNA, OverHangDNA]
    """
    return (
        OverHangDNA(fdna[: t_overhang.tprime], fdna[: -f_overhang.tprime]),
        OverHangDNA(tdna[: f_overhang.fprime], tdna[: -t_overhang.fprime]),
    )
