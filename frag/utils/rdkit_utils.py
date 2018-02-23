import os

from rdkit.Chem import ChemicalFeatures
from rdkit import Chem

# Generate ph4s for a molecule
class RDKitPh4(object):

    factory = None

    def __init__(self):
        # Generate the factory on init
        self.get_factory()

    def get_factory(self):
        """
        Generate the Ph4 feature factory
        :return:
        """
        if self.factory is None:
            this_dir, this_filename = os.path.split(__file__)
            data_path = os.path.join(this_dir, "data", "RDKitPh4.fdef")
            self.factory = ChemicalFeatures.BuildFeatureFactory(data_path)
        return self.factory


    def generate_ph4_for_mol(self, rdmol):
        """
        Generate a pharmacophore from an input molecule and a feature factory.
        :param rdmol: the input RDKit molecule
        :param factory: the feature factory
        :return: a list of 4 tuples (x,y,z, feature)
        """
        feats = self.get_factory().GetFeaturesForMol(rdmol)
        return [(feat.GetPos().x,feat.GetPos().y,feat.GetPos().z,feat.GetType()) for feat in feats]

def _parse_ligand_sdf(input_file):
    """
    Function to parse a series of ligands - return RDKit mols.
    :param input_file: the file to parse
    :param input_type: the type of ligands
    :return: the molecules parsed
    """
    mols = Chem.SDMolSupplier(input_file)
    return mols


def _get_c_of_mass(rdmol):
    """
    Get the unweighted centre of mass of an RDKit Molecule
    :param rdmol:
    :return:
    """
    atoms = rdmol.GetAtoms()
    conf = rdmol.GetConformer()
    x_coord = y_coord = z_coord = 0.0
    numatoms = 0.0
    # Assume all heavy atoms have the same mass
    for atom in atoms:
        if atom.GetAtomicNum() == 1 or atom.GetSmarts() == "[*]":
            continue
        numatoms += 1.0
        coords = conf.GetAtomPosition(atom.GetIdx())

        x_coord += float(coords.x)
        y_coord += float(coords.y)
        z_coord += float(coords.z)
    # Now we have all the coords -> we want to loop through
    if numatoms == 0:
        raise ValueError("No atoms in Molecules")
        return None, None, None
    return x_coord / numatoms, y_coord / numatoms, z_coord / numatoms


def _get_waters(lines):
    """Helper function to extract waters from a PDB file"""
    return [line for line in lines]# if line[17:20] == "HOH"]


def _get_water_coords(waters):
    """Helper function to get the coordinates from a load of waters."""
    rd_waters = Chem.MolFromPDBBlock("\n".join(waters))
    out_list = []
    if rd_waters is None:
        print("Warning - unable to parse waters.")
    if rd_waters is not None:
        # Check the waters exist
        conf = rd_waters.GetConformer()
        # Delete them for this protein
        for i in range(rd_waters.GetNumAtoms()):
            cp = conf.GetAtomPosition(i)
            if rd_waters.GetAtomWithIdx(i).GetSmarts() != "O":
                print("Warning - skipping a water")
                continue
            out_list.append((cp.x,cp.y,cp.z))
    return out_list


def _get_file(file_path,output_format,file_counter):
    if output_format=="smi":
        return Chem.SmilesWriter(file_path+"_"+str(file_counter)+".smi")
    else:
        return Chem.SDWriter(file_path+"_"+str(file_counter)+".sdf")

def _parse_mols(input_file,input_format):
    if input_format=="smi":
        return Chem.SmilesMolSupplier(input_file)
    else:
        return Chem.SDMolSupplier(input_file)

def _parse_pdb(data):
    return Chem.MolFromPDBFile(data)