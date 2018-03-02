from frag.alysis.cluster import dp_means
from frag.utils.parser import parse_ligand_ph4s

PH4_LAMBDA = 1.0
C_OF_M_LAMBDA = 6.0


def build_type_dict(mol_ph4_list,identifiers):
    type_dict = {}
    for i,mol in enumerate(mol_ph4_list):
        for ph4 in mol:
            x = ph4[0]
            y = ph4[1]
            z = ph4[2]
            ph4_type = ph4[3]
            if ph4_type in type_dict:
                type_dict[ph4_type]["coords"].append((x,y,z))
                type_dict[ph4_type]["mols"].append(identifiers[i])
            else:
                type_dict[ph4_type]={"coords":[(x,y,z)],"mols":[identifiers[i]]}
    return type_dict

def map_cluster(dp_means_cluster,mol_id_list):
    """

    :param dp_means_cluster:
    :param mol_id_list:
    :return:
    """
    out_dict = {}
    for cluster in dp_means_cluster.clusters:
        out_dict[cluster] = {"centre_of_mass": dp_means_cluster.clusters[cluster],
                             "mol_ids": []}
    for i,cluster_id in enumerate(dp_means_cluster.dataClusterId):
        out_dict[cluster_id]["mol_ids"].append(mol_id_list[i])
    return out_dict


def cluster_dp(vect_list, lam, mol_list):
    """
    Perform a DP Means clustering.
    :param vect_list: a list of lists of vectors
    :param lam: the clustering parameters
    :param mol_list: the molecular identifers in the same order as the list of vectors
    :return: a dictionary of the form {cluster_id: {centre_of_mass: (x,y,z), mol_ids: [1,5,12]}}
    """
    return map_cluster(dp_means(vect_list, lam), mol_list)


def run_lig_cluster(mols, identifiers):
    """
    Cluster a list of molecules
    :param mols: the input list of RDKit molecules
    :param identifiers: the corresponding list of identifiers
    :return:
    """
    # First we get the list of mols with their Ph4s
    mol_ph4_list = parse_ligand_ph4s(mols)
    # Then we build a dict of type: coords: [coords list], mols: [mol_index]
    type_dict = build_type_dict(mol_ph4_list,identifiers)
    # Then we cluster coords
    clusters = {}
    for ph4_type in type_dict:
        if ph4_type == "c_of_m":
            clusters[ph4_type]= cluster_dp(type_dict[ph4_type]["coords"],C_OF_M_LAMBDA,type_dict[ph4_type]["mols"])
        else:
            clusters[ph4_type]= cluster_dp(type_dict[ph4_type]["coords"],PH4_LAMBDA,type_dict[ph4_type]["mols"])
    return clusters