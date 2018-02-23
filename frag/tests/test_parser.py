import unittest
from frag.utils import parser



class ParserTest(unittest.TestCase):
    water_data = """HETATM 2008  O   HOH B 184      53.034 -39.489  96.872  1.00 67.70           O
HETATM 2010  O   HOH B 186      39.366 -30.950  88.735  1.00 66.27           O
HETATM 2011  O   HOH B 187      38.861 -67.134  82.852  1.00 69.11           O
HETATM 2012  O   HOH B 188      48.438 -40.466  97.529  1.00 41.98           O
HETATM 2015  O   HOH B 190      47.858 -60.571  77.866  1.00 52.55           O
HETATM 2016  O   HOH B 191      52.415 -50.993  68.148  1.00 50.73           O
HETATM 2017  O   HOH B 192      48.922 -42.540  98.150  1.00 52.43           O
HETATM 2018  O   HOH B 193      60.968 -55.453  92.185  1.00 37.56           O
HETATM 2019  O   HOH B 194      26.058 -55.837  68.104  1.00 60.13           O
HETATM 2020  O   HOH B 195      26.923 -52.747  83.917  1.00 59.04           O
HETATM 2021  O   HOH B 196      42.376 -40.566  71.629  1.00 45.71           O
HETATM 2022  O   HOH B 197      45.966 -45.361 100.995  1.00 56.39           O
HETATM 2023  O   HOH B 198      40.498 -49.338  59.760  1.00 74.54           O
HETATM 2024  O   HOH B 199      62.006 -56.842  90.642  1.00 50.69           O"""
    single_water = """HETATM 2008  O   HOH B 184      53.034 -39.489  96.872  1.00 67.70           O"""

    def test_water_parser(self):
        out_data = parser._get_waters(self.water_data.split("\n"))
        self.assertEqual(len(out_data),14)
        out_data = parser._get_waters(self.single_water.split("\n"))
        self.assertEqual(len(out_data),1)

    def test_water_reader(self):
        out_data = parser._get_waters(self.water_data.split("\n"))
        water_coords = parser._get_water_coords(out_data)
        self.assertEqual(len(water_coords),14)
        self.assertAlmostEqual(water_coords[4][2],77.866)
        out_data = parser._get_waters(self.single_water.split("\n"))
        water_coords = parser._get_water_coords(out_data)
        self.assertEqual(len(water_coords), 1)
        self.assertAlmostEqual(water_coords[0][1],-39.489)


if __name__ == '__main__':
    unittest.main()