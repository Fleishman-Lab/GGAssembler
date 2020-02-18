import fire
from dawdlib.golden_gate.find_gg import create_goldengates as cg_main
from dawdlib.create_embl.create_embl import create_embl


def create_goldengates(
    dna_file: str,
    deg_table_file: str,
    out_dir: str,
    min_var_oligo_length: int,
    max_var_oligo_length: int,
    min_const_oligo_length: int,
    no_solutions: int = 1,
    min_oligos: int = None,
    max_oligos: int = None,
    gate_self_binding_min: int = 2000,
    gate_crosstalk_max: int = 1000,
):
    cg_main(
        dna_file,
        deg_table_file,
        out_dir,
        min_var_oligo_length,
        max_var_oligo_length,
        min_const_oligo_length,
        no_solutions,
        min_oligos,
        max_oligos,
        gate_self_binding_min,
        gate_crosstalk_max,
    )


def embl(deg_table, embl_in, gate_path, embl_out: str):
    create_embl(
        deg_table_file=deg_table,
        embl_file=embl_in,
        path_file=gate_path,
        out_embl_file=embl_out,
    )


if __name__ == "__main__":
    fire.Fire()
