import pandas as pd


def load_deg_table(deg_table_path: str) -> pd.DataFrame:
    return pd.read_csv(deg_table_path, na_filter=True, keep_default_na=False,)


def parse_degenerate_codon_csv(csv_file: str) -> pd.DataFrame:
    """
    parse the degenerate codons table
    Args:
        csv_file (str): the csv file for the degenerate codons

    Returns:
        (pd.DataFrame) describing the required degenerate codons

    """
    return pd.read_csv(
        csv_file,
        index_col=False,
        na_values="NaN",
        converters={
            "ENCODED_AAS": lambda x: x.strip("[]").replace("'", "").split(", "),
            "ENCODED_COUNT": lambda x: [
                int(a) for a in x.strip("[]").replace("'", "").split(", ")
            ],
        },
    )
