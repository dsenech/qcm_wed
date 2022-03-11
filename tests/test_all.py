import unittest
import os

DIR_PATH = os.path.dirname(os.path.realpath(__file__))


def run_file(test_name):
    """Runs a file of a given name inside the testing directory. Mostly aesthetic.
    """
    return os.system(f"python {DIR_PATH}/{test_name}")


class TestAll(unittest.TestCase):

    def test_averages_bath(self):
        self.assertFalse(run_file("test_averages_bath.py"))

    def test_averages(self):
        self.assertFalse(run_file("test_averages.py"))

    def test_berry(self):
        self.assertFalse(run_file("test_berry.py"))

    def test_cdmft(self):
        self.assertFalse(run_file("test_cdmft.py"))

    def test_fixed_density_loop(self):
        self.assertFalse(run_file("test_fixed_density_loop.py")) ######### needs to be shortened

    def test_hybridization(self):
        self.assertFalse(run_file("test_hybridization.py"))

    def test_instances(self):
        self.assertFalse(run_file("test_instances.py"))

    def test_mixing_cf(self):
        self.assertFalse(run_file("test_mixing_cf.py"))

    def test_mixing(self):
        self.assertFalse(run_file("test_mixing.py"))

    def test_mixing2(self):
        self.assertFalse(run_file("test_mixing2.py"))

    def test_spectral(self):
        self.assertFalse(run_file("test_spectral.py")) ################# make sure that all is output as pdf to allow user-less execution (applies to the previous mixing I think)

    def test_vca(self):
        self.assertFalse(run_file("test_vca.py")) ####### graphs need to go to pdf


if __name__ == "__main__":
    unittest.main()