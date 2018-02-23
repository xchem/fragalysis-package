from frag.alysis.cluster import dp_means
from frag.utils.parser import parse_waters, parse_ligand_ph4s, parse_ligands, parse_residues
from frag.alysis.models import ClusterStuff

CLUSTER_DICT = {"waters": ClusterStuff(parser=parse_waters, lamb=0.7, cluster=dp_means),
                "ph4": ClusterStuff(parser=parse_ligand_ph4s, lamb=0.7, cluster=dp_means),
                "residues": ClusterStuff(parser=parse_residues, lamb=0.7, cluster=dp_means),
                "sites": ClusterStuff(parser=parse_ligands, lamb=0.7, cluster=dp_means)}


def run_cluster(input_list):
    # Residues
    CLUSTER_DICT["residues"].run([(x.resid_pdb,x.struct_id) for x in input_list])
    # Waters
    CLUSTER_DICT["waters"].run([(x.water_pdb,x.struct_id) for x in input_list])
    # Pharmacophores
    CLUSTER_DICT["ph4"].run([(x.ligand,x.struct_id) for x in input_list])
    # Sites
    CLUSTER_DICT["sites"].run([(x.ligand,x.struct_id) for x in input_list])