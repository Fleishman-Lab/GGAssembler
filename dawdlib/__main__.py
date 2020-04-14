import fire

from dawdlib.golden_gate.reaction_sim import ReactionCLI


class CLI:
    def __init__(self):
        self.ggsim = ReactionCLI()

    @staticmethod
    def degenerate_table(
        res_file: str, output_csv: str, organism_id: str = "37762"
    ) -> None:
        """
        Generate the degenerate codons required in order to create the diversity specified in the res file.
        Args:
            res_file (str): Input res file to find degenerate AA codons.
            output_csv (str): Output file location.
            organism_id (str): Organism ID, default: Escherichia Coli. Used to find the best codon for each AA.
        """
        from dawdlib.degenerate_dna.deg_table import generate_deg_csv

        generate_deg_csv(res_file, output_csv, organism_id)

    @staticmethod
    def create_goldengates(
        dna_file: str,
        deg_table_file: str,
        out_dir: str,
        min_var_oligo_length: int,
        max_var_oligo_length: int,
        min_const_oligo_length: int,
        gg_temp: int,
        gg_time: int,
        gg_min_efficiency: float,
        gg_min_fidelity: float,
        no_solutions: int = 1,
        min_oligos: int = None,
        max_oligos: int = None,
    ) -> None:
        """

        Args:
            dna_file (str):
            deg_table_file (str):
            out_dir (str):
            min_var_oligo_length (int):
            max_var_oligo_length (int):
            min_const_oligo_length (int):
            gg_temp:
            gg_time:
            gg_min_efficiency:
            gg_min_fidelity:
            no_solutions (int):
            min_oligos (int):
            max_oligos (int):

        Returns:
            object:
        """
        from dawdlib.golden_gate.find_gg import create_goldengates as cg_main

        cg_main(
            dna_file,
            deg_table_file,
            out_dir,
            min_var_oligo_length,
            max_var_oligo_length,
            min_const_oligo_length,
            gg_temp,
            gg_time,
            gg_min_efficiency,
            gg_min_fidelity,
            no_solutions,
            min_oligos,
            max_oligos,
        )

    @staticmethod
    def embl(embl_in: str, embl_out: str, deg_table: str = "", gate_path: str = ""):
        """

        ArgsReactionCLI:
            embl_in (str): file in embl format to write features to
            embl_out (str): file where to write embl file with all features
            deg_table (str): (optional) degenerate codon table file in csv format
            gate_path (str): (optional) file with list of gates
        """
        from dawdlib.embl_utils.create_embl import create_embl

        create_embl(
            deg_table_file=deg_table,
            embl_file=embl_in,
            path_file=gate_path,
            out_embl_file=embl_out,
        )

    @staticmethod
    def combine(
        dc_table: str,
        gate_path: str,
        dna: str,
        to_order: str,
        prefix: str = "CGTGCGGTCTCG",
        suffix: str = "CGAGACCGCGCCGGGC",
    ):
        """
        combines the degenerate codons table and the path of gates to create the oligos required to create the library.

        Args:
            dc_table: degenerate codon table file in csv format
            gate_path (str): file with list of gates
            dna (str): file with the full gene DNA
            to_order (str): file path to write the table of oligos to order
            prefix (str): the string to add before each oligo. defaults to a BsaI recognition site
            suffix (str): the string to add after each oligo. defaults to a BsaI recognition site and primer
        """
        from dawdlib.gg_dc_combine.gg_dc_combine import combine_gate_path_deg_codons

        combine_gate_path_deg_codons(
            dc_table_file=dc_table,
            gate_path_file=gate_path,
            dna_file=dna,
            prefix=prefix,
            suffix=suffix,
            to_order_df_file=to_order,
        )


def main(args=None):
    fire.Fire(CLI, name="dawdlib")


if __name__ == "__main__":
    main()
