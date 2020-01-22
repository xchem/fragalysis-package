import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit import RDConfig


fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
        'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')


def get_FeatureMapScore(small_m, large_m, score_mode=FeatMaps.FeatMapScoreMode.All):
    featLists = []
    for m in [small_m, large_m]:
        rawFeats = fdef.GetFeaturesForMol(m)
        # filter that list down to only include the ones we're intereted in
        featLists.append([f for f in rawFeats if f.GetFamily() in keep])
    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
    fms[0].scoreMode = score_mode
    fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
    return fm_score


def score(reflig, prb_mols, ids, score_mode=FeatMaps.FeatMapScoreMode.All, p=False):
    ref = Chem.AddHs(reflig)
    idx = 0

    results_sucos = {}
    results_tani = {}

    smi_mol = Chem.MolToSmiles(prb_mols)

    for i in ids:

        prb = Chem.AddHs(Chem.MolFromMolBlock(Chem.MolToMolBlock(prb_mols, confId=i)))

        fm_score = get_FeatureMapScore(ref, prb, score_mode)
        fm_score = np.clip(fm_score, 0, 1)

        protrude_dist = rdShapeHelpers.ShapeProtrudeDist(ref, prb,
                                                         allowReordering=False)
        protrude_dist = np.clip(protrude_dist, 0, 1)

        SuCOS_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)
        tanimoto_score = Chem.rdShapeHelpers.ShapeTanimotoDist(ref, prb)

        results_sucos[str(idx)] = SuCOS_score
        results_tani[str(idx)] = tanimoto_score

        if p:
            print("********************************")
            print("index: " + str(idx))
            print("SuCOS score:\t%f" % SuCOS_score)
            print("Tani score:\t%f" % tanimoto_score)
            print("********************************")

        idx += 1

    return results_sucos


def get_best_align(hit_mblock, elab_smiles):
    hit_mol = Chem.MolFromMolBlock(hit_mblock)
    elab_mol = Chem.MolFromSmiles(elab_smiles)
    ids = AllChem.EmbedMultipleConfs(elab_mol, numConfs=100, params=AllChem.ETKDG())

    for cid in ids:
        o3d = Chem.rdMolAlign.GetO3A(prbMol=elab_mol, refMol=hit_mol, prbCid=cid)
        o3d.Align()

    results_sucos = score(hit_mol, elab_mol, ids)
    best_i = list(results_sucos.values()).index(max(results_sucos.values()))
    elab_molblock = Chem.MolToMolBlock(elab_mol, confId=best_i)

    return elab_molblock


def get_core_mol(three_d_mol, core_mol):
    repl_sidechains = AllChem.ReplaceSidechains(three_d_mol, core_mol)
    if repl_sidechains:
        return AllChem.DeleteSubstructs(repl_sidechains, Chem.MolFromSmiles("*"))


def gen_conf_from_vector(input_mol_block, elaborated_smiles):
    # Get the mol
    m = get_best_align(input_mol_block, elaborated_smiles)
    return m

def generate_confs_for_vector(input_smiles, input_mol_block):
    confs = []
    for elaborated_smiles in input_smiles:
        confs.append(
            gen_conf_from_vector(input_mol_block, elaborated_smiles)
        )
    return confs


def generate_confs_for_vectors(input_vectors, input_mol_block):
    confs = {}
    for vector in input_vectors:
        confs[vector] = generate_confs_for_vector(
            input_vectors[vector], input_mol_block
        )
    return confs


if __name__ == "__main__":
    conf = generate_confs_for_vectors(
        {"O=C(Nc1ccc([Xe])nc1)N1CCNCC1_1_ADDITION": ["COc1ccc(NC(=O)N2CCNCC2)cn1"]},
        "\n     RDKit          3D\n\n 15 16  0  0  0  0  0  0  0  0999 V2000\n   -8.0100   -6.1770   18.7960 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -7.8350   -7.3610   19.1060 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.7970   -8.0980   18.5710 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.0500   -7.8150   17.3890 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.2220   -6.6910   16.5800 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.5070   -6.4830   15.4790 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.5740   -7.3730   15.1470 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.3180   -8.4870   15.8900 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.0700   -8.7240   17.0200 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -8.6530   -7.9840   20.0040 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -8.8660   -9.4290   20.1410 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -10.3240   -9.8000   19.9120 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -11.2290   -8.6780   20.2340 N   0  0  0  0  0  0  0  0  0  0  0  0\n  -10.7100   -7.8960   21.3670 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -9.4210   -7.2020   20.9820 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  2 10  1  0\n  3  4  1  0\n  4  5  2  0\n  4  9  1  0\n  5  6  1  0\n  6  7  2  0\n  7  8  1  0\n  8  9  2  0\n 10 11  1  0\n 10 15  1  0\n 11 12  1  0\n 12 13  1  0\n 13 14  1  0\n 14 15  1  0\nM  END\n",
    )

    print(conf)
