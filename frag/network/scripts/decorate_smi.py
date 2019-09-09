from frag.network.decorate import decorate_smi
import argparse
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from tqdm import tqdm
from frag.network.models import Attr

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Decorate a library of molecules for insertion to the database.')
    parser.add_argument('--input_smi')
    parser.add_argument('--output_attr')
    args = parser.parse_args()
    out_smi = open(args.output_attr,"w")
    for mol in tqdm(Chem.SmilesMolSupplier(args.input_smi,delimiter=',',smilesColumn=1,nameColumn=0)):
        this_smi = Chem.MolToSmiles(mol,isomericSmiles=True)
        new_smis = decorate_smi(this_smi)
        new_murck = decorate_smi(MurckoScaffold.MurckoScaffoldSmiles(this_smi))
        # mol_frags = get_fragments(Chem.MolFromSmiles(this_smi),iso_labels=False)
        # new_smis.extend([x.replace("Xe","At") for x in mol_frags])
        new_smis.extend(new_murck)
        new_smis = list(set(new_smis))
        # Do this on original and on Murcko Scaffold
        name = mol.GetProp("_Name")
        new_attr = Attr(this_smi, ["EM", name])
        out_smi.write(str(new_attr) + "\n")
        for i,smi in enumerate(new_smis):
            new_attr = Attr(smi,["EM",name+"_"+str(i)])
            out_smi.write(str(new_attr)+"\n")