from rdkit import Chem
from rdkit.Chem import AllChem


def get_core_mol(three_d_mol, core_mol):
    repl_sidechains = AllChem.ReplaceSidechains(three_d_mol, core_mol)
    if repl_sidechains:
        return AllChem.DeleteSubstructs(repl_sidechains, Chem.MolFromSmiles("*"))


def gen_conf_from_vector(input_mol_block, vector, elaborated_smiles):
    # Get the mol
    output_mol = Chem.MolFromSmiles(elaborated_smiles)
    three_d_mol = Chem.MolFromMolBlock(input_mol_block)
    # Get the shared core
    core = get_core_mol(
        three_d_mol,
        Chem.RemoveHs(Chem.MolFromSmiles(vector.split("_")[0].replace("Xe", "H"))),
    )
    if core:
        return Chem.MolToMolBlock(AllChem.ConstrainedEmbed(output_mol, core))


def generate_confs_for_vector(input_vector, input_smiles, input_mol_block):
    confs = []
    for elaborated_smiles in input_smiles:
        confs.append(
            gen_conf_from_vector(input_mol_block, input_vector, elaborated_smiles)
        )
    return confs


def generate_confs_for_vectors(input_vectors, input_mol_block):
    confs = {}
    for vector in input_vectors:
        confs[vector] = generate_confs_for_vector(
            vector, input_vectors[vector], input_mol_block
        )
    return confs


if __name__ == "__main__":
    generate_confs_for_vectors(
        {"O=C(Nc1ccc([Xe])nc1)N1CCNCC1_1_ADDITION": ["COc1ccc(NC(=O)N2CCNCC2)cn1"]},
        "\n     RDKit          3D\n\n 15 16  0  0  0  0  0  0  0  0999 V2000\n   -8.0100   -6.1770   18.7960 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -7.8350   -7.3610   19.1060 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.7970   -8.0980   18.5710 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.0500   -7.8150   17.3890 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.2220   -6.6910   16.5800 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.5070   -6.4830   15.4790 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.5740   -7.3730   15.1470 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.3180   -8.4870   15.8900 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.0700   -8.7240   17.0200 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -8.6530   -7.9840   20.0040 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -8.8660   -9.4290   20.1410 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -10.3240   -9.8000   19.9120 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -11.2290   -8.6780   20.2340 N   0  0  0  0  0  0  0  0  0  0  0  0\n  -10.7100   -7.8960   21.3670 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -9.4210   -7.2020   20.9820 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  2 10  1  0\n  3  4  1  0\n  4  5  2  0\n  4  9  1  0\n  5  6  1  0\n  6  7  2  0\n  7  8  1  0\n  8  9  2  0\n 10 11  1  0\n 10 15  1  0\n 11 12  1  0\n 12 13  1  0\n 13 14  1  0\n 14 15  1  0\nM  END\n",
    )
