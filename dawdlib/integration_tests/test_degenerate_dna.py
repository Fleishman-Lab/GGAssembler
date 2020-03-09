import os
import unittest


class IntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.here = os.path.dirname(__file__)
        self.deg_table_file = os.path.join(self.here, "degenerate_table.csv")
        self.resfile_file = os.path.join(
            self.here, "resources", "271_2p48_reduced.resfile"
        )
        self.dna_file = os.path.join(self.here, "resources", "271_2p48.fasta")

    def test_degenerate_table(self):
        print(self.resfile_file)
        from dawdlib.degenerate_dna.deg_table import generate_deg_csv

        generate_deg_csv(self.resfile_file, self.deg_table_file, "37762")

    def test_create_goldengates(self):
        pass


if __name__ == "__main__":
    unittest.main()
