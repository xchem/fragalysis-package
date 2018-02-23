from frag.utils.rdkit_utils import RDKitPh4,_get_c_of_mass,find_dist,_get_res_rmsds,get_res_atom_name,_get_res
import unittest
from rdkit import Chem

SDF_DATA = """
     RDKit          3D

  8  8  0  0  0  0  0  0  0  0999 V2000
   -0.2375    1.2479   -0.2221 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5033    1.0878    0.3249 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9231   -0.1960    0.6598 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4265   -1.1855   -0.1809 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0741   -1.1112   -0.4792 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5447    0.1297   -0.4205 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0371    0.2224   -0.4859 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5826   -0.1950    0.8038 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  1  0
  6  1  1  0
M  END

"""

PDB_DATA = """ATOM   1315  N   HIS A 174      29.516 -40.779  92.218  1.00 59.90           N
ATOM   1316  CA  HIS A 174      30.777 -40.579  92.954  1.00 55.77           C
ATOM   1317  C   HIS A 174      31.928 -40.122  92.037  1.00 48.90           C
ATOM   1318  O   HIS A 174      32.167 -40.710  90.977  1.00 53.82           O
ATOM   1319  CB  HIS A 174      31.169 -41.858  93.717  1.00 59.44           C
ATOM   1320  CG  HIS A 174      30.133 -42.338  94.691  1.00 63.61           C
ATOM   1321  ND1 HIS A 174      29.504 -41.501  95.589  1.00 63.97           N
ATOM   1322  CD2 HIS A 174      29.645 -43.580  94.931  1.00 66.30           C
ATOM   1323  CE1 HIS A 174      28.658 -42.202  96.324  1.00 66.58           C
ATOM   1324  NE2 HIS A 174      28.725 -43.467  95.945  1.00 69.74           N
"""

class Ph4Test(unittest.TestCase):

    def test_ph4_parser(self):
        """
        Test the pharmacophore generation.
        Need more test cases
        :return:
        """
        rdmol = Chem.MolFromMolBlock(SDF_DATA)
        rdkit_ph4 = RDKitPh4()
        feats = rdkit_ph4.generate_ph4_for_mol(rdmol=rdmol)
        self.assertEqual(len(feats),5)
        self.assertLessEqual([x[3] for x in feats],
                             ['SingleAtomDonor', 'BasicGroup', 'Arom6', 'ThreeWayAttach', 'RH6_6'])
        self.assertAlmostEqual([x[1] for x in feats][1],
                               -0.195)

class CentreOfMassTest(unittest.TestCase):
    def test_centre_of_mass(self):
        rdmol = Chem.MolFromMolBlock(SDF_DATA)
        centre_of_mass = _get_c_of_mass(rdmol)
        self.assertAlmostEqual(centre_of_mass[0],-1.2499999999970868e-05)
        self.assertAlmostEqual(centre_of_mass[1],1.2499999999995154e-05)
        self.assertAlmostEqual(centre_of_mass[2],-1.2499999999984746e-05)

class MathTest(unittest.TestCase):

    def test_square_dist(self):
        self.assertEqual(find_dist(8,22,12,4,6,6),308.0)

class ResTest(unittest.TestCase):

    def test_get_res_rmsd(self):
        input_data = [
                       {"A": [1.0, 2.0, 3.0],
                        "B": [1.0, 2.0, 3.0],
                        "C": [1.0, 2.0, 3.0],
                        },
                       {"A": [2.0, 3.0, 4.0],
                        "B": [1.0, 2.0, 3.0],
                        "C": [4.0, 3.0, 2.0],
                        },
                        {"A": [-0.2375, 1.2479, -0.2221],
                         "B": [-1.5033, 1.0878, 0.3249],
                         "C": [-1.9231, -0.1960, 0.6598],
                        }
        ]
        output_data = [[2.160246899469287, 3.8977444101257674],
                       [2.160246899469287, 5.392922351255084],
                       [3.8977444101257674, 5.392922351255084]]
        test_output = _get_res_rmsds(input_data)
        print(test_output)
        for i,output in enumerate(test_output):
            self.assertEqual(len(output),2)
            self.assertEqual(output_data[i][0],output[0])

    def test_get_res_name(self):
        pdb_mol = Chem.MolFromPDBBlock(PDB_DATA)
        conf = pdb_mol.GetConformer()
        atoms = pdb_mol.GetAtoms()
        for atom in atoms:
            atom
        unique_name, atom_name, position = get_res_atom_name(atom,conf)
        self.assertEqual(unique_name,"174_A_HIS")
        self.assertEqual(atom_name,"NE2")
        self.assertAlmostEqual(position[0],28.725)
        self.assertAlmostEqual(position[1],-43.467)
        self.assertAlmostEqual(position[2],95.945)

    def test_get_res(self):
        mol  = Chem.MolFromPDBBlock(PDB_DATA)
        out_dict = _get_res(mol)
        self.assertListEqual(out_dict.keys(),["174_A_HIS"])