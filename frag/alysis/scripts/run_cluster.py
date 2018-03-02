from frag.alysis.run_clustering import run_cluster
from frag.alysis.models import StructHolder
### Example of how it would work

def run_cluster_proasis():

    # This can be collated by querying RestAPI. Actual input would be list of struct_ids - or target.
    proasis_input = [StructHolder(struct_id="struct_one",apo_pdb="struct_one_apo.pdb", ligand="struct_one_ligand.sdf"),
                     StructHolder(struct_id="struct_two", apo_pdb="struct_two_apo.pdb", ligand="struct_two_ligand.sdf")]
    run_cluster(input_list=proasis_input)

# This would be input ("struct_one_apo.pdb","struct_one_ligand.sdf","struct_one")
### Need script that splits waters from residues (based on group). Possbily parsing with PDB


def run_cluster_simple():
    ## Create this test data
    input_list = [StructHolder(struct_id="struct_one",resid_pdb="struct_one_resid.pdb", ligand="struct_one_ligand.sdf",
                               water_pdb="struct_one_water.pdb"),
                  StructHolder(struct_id="struct_two", resid_pdb="struct_two_resid.pdb", ligand="struct_two_ligand.sdf",
                               water_pdb="struct_two_water.pdb"),]
    run_cluster(input_list=input_list)
