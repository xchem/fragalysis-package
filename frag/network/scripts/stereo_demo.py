# Example of standardization and iso/noniso smiles
# This script is only for demonstration purposes.


from rdkit import Chem
from frag.utils.rdkit_utils import standardize


smiles1='N[C@@H](CC1=CC=CC=C1)C(O)=O'            # an isomeric molecule
smiles2='CC1=CC=CC=C1'                           # a non-isomeric molecule
smiles3 = '[Na+].N[C@@H](CC1=CC=CC=C1)C([O-])=O' # molecule needing standardization

mol1 = Chem.MolFromSmiles(smiles1)
mol2 = Chem.MolFromSmiles(smiles2)
mol3 = Chem.MolFromSmiles(smiles3)


print("\nIsomeric molecule 1")
iso1 = Chem.MolToSmiles(mol1, isomericSmiles=True, canonical=True)
noniso1 = Chem.MolToSmiles(mol1, isomericSmiles=False, canonical=True)

print("Read:        " + smiles1)
print("Isomeric:    " + iso1)
print("NonIsomeric: " + noniso1)
print("Is isomeric? " + str(iso1 != noniso1)) # iso and noniso smiles are different


print("\nIsomeric molecule 2")
iso2 = Chem.MolToSmiles(mol2, isomericSmiles=True, canonical=True)
noniso2 = Chem.MolToSmiles(mol2, isomericSmiles=False, canonical=True)

print("Read:        " + smiles2)
print("Isomeric:    " + iso2)
print("NonIsomeric: " + noniso2)
print("Is isomeric? " + str(iso2 != noniso2)) # iso and noniso smiles are the same


print("\nNon-standardized molecule 3")
std3 = standardize(mol3) # gives back a RDKit RWMol
iso3 = Chem.MolToSmiles(std3, isomericSmiles=True, canonical=True)
noniso3 = Chem.MolToSmiles(std3, isomericSmiles=False, canonical=True)

print("Read:        " + smiles3) # standardization changes the smiles
print("Isomeric:    " + iso3)
print("NonIsomeric: " + noniso3)
print("Is isomeric? " + str(iso3 != noniso3))
print("Is smiles1?  " + str(iso1 == iso3))