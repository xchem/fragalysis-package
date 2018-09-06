from neomodel import (
    StructuredNode,
    StructuredEdge,
    StringProperty,
    FloatProperty,
    IntegerProperty,
    RelationshipTo,
    RelationshipFrom,
    BooleanProperty,
)


class F2Edge(StructuredEdge):
    exclude_smiles = StringProperty()
    exclude_type = StringProperty()
    rebuilt_smiles = StringProperty()
    rebuilt_ring_smiles = StringProperty()
    rebuilt_type = StringProperty()
    excluded_ring_smiles = StringProperty()
    # The Exclude and Rebuilt smiles combined
    both_sides = StringProperty(unique_index=True)


class BaseMolecule(StructuredNode):
    smiles = StringProperty(unique_index=True)
    heavy_atom_count = IntegerProperty()
    ring_atom_count = IntegerProperty()
    ring_smiles = StringProperty()
    # Annotations of price and molecule id
    price = IntegerProperty()
    mol_id = StringProperty()


class Compound(BaseMolecule):
    """
    A molecule - a real purchasable molecule
    """
    fragments = RelationshipTo("Fragment", "F2EDGE")
    conformers = RelationshipTo("Conformer", "CONFEDGE")
    parent_molecules = RelationshipTo("Compound", "PARENTEDGE")


class Fragment(BaseMolecule):
    """
    A fragment - not a real molecule
    """
    molecules = RelationshipFrom(Compound, "F2EDGE", model=F2Edge)


class Conformer(StructuredNode):
    """
    The 3D conformer of a molecule - either Docked or imaginary
    """
    # Going to go with uuid - seems as good as anything - this connects to the Django database for 3D information
    unique_string = StringProperty(unique_index=True)
    conformer = RelationshipFrom(Compound, "CONFEDGE")
    score = FloatProperty()
    is_xtal_pose = BooleanProperty()
    threed_sim_to_reference = FloatProperty()
    rmsd_to_xtal_pose = FloatProperty()
    # Connect to Molecule and Protein


class ConformerFeature(StructuredNode):
    """
    A pharmacophore feature for a given Conformer
    """
    # Unique id - uuid:PHARMAFEATURE:INDEX
    unique_string = StringProperty(unique_index=True)
    pharma_feature = StringProperty()
    index = IntegerProperty()
    # Connect to SpecificResidue and to other SpecificResidueFeature annotated by distance and angle


class Target(StructuredNode):
    """
    A target is defined as a unique Uniprot accessiion
    """
    # A given uniprot accession
    uniprot_id = StringProperty(unique_index=True)
    # Now a series of properties
    sequence = StringProperty()
    # Connect Target-Target (sequence similarity data from Biojava)


class Protein(StructuredNode):
    """
    A unique set of 3D coordinates for a protein structure
    """
    # A given protein structure
    sequence = StringProperty()
    # Connect Protein-Target and Protein-Protein annotate by they're sequence identity


class GeneralResidue(StructuredNode):
    """
    A unique residue type - as defined by PDB RES codes
    """
    # A particular type of residue (defined by PDB RES code)
    res_code = StringProperty(unique_index=True)
    # Connect Residue to Residue based on their similarity. Human defined at first


class SpecificResidue(StructuredNode):
    """
    A specific 3D coords for a given residue
    """
    # Unique id - PDBID:CHAINID:RES_NUM:ALTLOC
    unique_string = StringProperty(unique_index=True)
    pdb_id = StringProperty()
    chain_id = StringProperty()
    residue_number = StringProperty()
    res_id = StringProperty()
    # Connect SpecificResidue to GeneralResidue and Protein
    # and SpecificResidue to SpecificResidue based on orientation and distance


class SpecificResidueFeature(StructuredNode):
    """
    A pharmacophore feature for a given SpecificResidue
    """
    # Unique id - PDBID:CHAINID:RES_NUM:ALTLOC:PHARMATYPE:INDEX
    unique_string = StringProperty(unique_index=True)
    # The 3 Letter code residue id
    res_id = StringProperty()
    pharma_type = StringProperty()
    index = IntegerProperty()
    # Connect to SpecificResidue and
    # to other SpecificResidueFeature/ConformerFeature annotated by distance and angle
