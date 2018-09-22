import random
from frag.utils.network_utils import write_results, get_driver, canon_input


class ReturnObject(object):

    def __init__(self, start_smi, end_smi, label, edge_count, change_frag, iso_label):
        """
        Build this object.
        :param start_smi:
        :param end_smi:
        :param label:
        :param frag_type:
        :param edge_count:
        """
        self.start_smi = start_smi
        self.end_smi = end_smi
        self.label = label
        self.iso_label = iso_label
        self.frag_type = None
        self.edge_count = edge_count
        self.change_frag = change_frag

    def __str__(self):
        out_list = [self.label, str(self.edge_count), self.frag_type]
        return "_".join(out_list)


def find_double_edge(tx, input_str):
    return tx.run(
        "MATCH (sta:F2 {smiles:$smiles})-[nm:F2EDGE]-(mid:F2)-[ne:F2EDGE]-(end:EM) where"
        " abs(sta.hac-end.hac) <= 3 and abs(sta.chac-end.chac) <= 1"
        " and sta.smiles <> end.smiles "
        "RETURN sta, nm, mid, ne, end "
        "order by split(nm.label, '|')[4], split(ne.label, '|')[2];",
        smiles=input_str,
    )


def find_triple_edge_growth(
    tx,
    input_str,
    heavy_atom_diff_min=6,
    heavy_atom_diff_max=10,
    mid_heavy_atom_diff_min=-1,
    mid_heavy_atom_diff_max=3,
):
    return tx.run(
        "MATCH (sta:F2 {smiles:$smiles})-[nm:F2EDGE]-(mid_one:F2)-[ne:F2EDGE]-(mid:EM)-[nm2:F2EDGE]-(mid_two:F2)-[ne2:F2EDGE]-(end:EM) where"
        " end.hac-sta.hac > $hacmin and end.hac-sta.hac <= $hacmax"
        " and mid.hac-sta.hac > $chacmin and mid.hac-sta.hac <= $chacmax"
        " and sta.smiles <> mid.smiles and sta.smiles <> end.smiles "
        " WITH collect("
        "{"
        "end: end.smiles,"
        "mid: mid.smiles,"
        "frag_one: mid_one.smiles,"
        "frag_two: mid_two.smiles"
        "}"
        ") AS edges"
        " RETURN edges",
        smiles=input_str,
        hacmin=heavy_atom_diff_min,
        hacmax=heavy_atom_diff_max,
        chacmin=mid_heavy_atom_diff_min,
        chacmax=mid_heavy_atom_diff_max,
    )


def add_follow_ups(tx, input_str):
    return tx.run(
        "MATCH (sta:F2 {smiles:$smiles})-[nm:F2EDGE]-(mid:F2)-[ne:F2EDGE]-(end:EM) where"
        " abs(sta.hac-end.hac) <= 3 and abs(sta.chac-end.chac) <= 1"
        " and sta.smiles <> end.smiles "
        " MERGE (end)-[:FOLLOW_UP]->(sta)",
        smiles=input_str,
    )


def find_proximal(tx, input_str):
    return tx.run(
        "match p = (n:F2{smiles:$smiles})-[nm]-(m:EM)"
        "where abs(n.hac-m.hac) <= 3 and abs(n.chac-m.chac) <= 1 "
        "return n, nm, m "
        "order by split(nm.label, '|')[4];",
        smiles=input_str,
    )


def find_custom(tx, input_str):
    return tx.run(input_str)


def get_type(r_group_form, sub_one, sub_two):
    if "." in r_group_form:
        if "C1" in sub_two:
            return "ring_linker"
        return "linker"
    if "C1" in sub_two:
        return "ring_replacement"
    return "replacement"


def define_double_edge_type(record):
    mol_one = record["sta"]
    label = str(record["ne"]["label"].split("|")[4])
    iso_label = str(record["ne"]["label"].split("|")[5])
    change_frag = str(record["ne"]["label"].split("|")[2])
    mol_two = record["mid"]
    mol_three = record["end"]
    diff_one = mol_one["hac"] - mol_two["hac"]
    diff_two = mol_two["hac"] - mol_three["hac"]
    ret_obj = ReturnObject(
        mol_one["smiles"], mol_three["smiles"], label, 2, change_frag, iso_label
    )
    if "." in label:
        ret_obj.frag_type = "LINKER"
    elif diff_one >= 0 and diff_two >= 0:
        ret_obj.frag_type = "DELETION"
    elif diff_one <= 0 and diff_two <= 0:
        ret_obj.frag_type = "ADDITION"
    else:
        ret_obj.frag_type = "REPLACE"
    return ret_obj


def define_proximal_type(record):
    """
    Define the type returned for proximal systems
    :param record:
    :return:
    """
    mol_one = record["n"]
    label = str(record["nm"]["label"].split("|")[4])
    iso_label = str(record["nm"]["label"].split("|")[5])
    change_frag = str(record["nm"]["label"].split("|")[2])
    mol_two = record["m"]
    ret_obj = ReturnObject(
        mol_one["smiles"], mol_two["smiles"], label, 1, change_frag, iso_label
    )
    if "." in label:
        ret_obj.frag_type = "LINKER"
    elif mol_one["hac"] - mol_two["hac"] > 0:
        ret_obj.frag_type = "DELETION"
    elif mol_one["hac"] - mol_two["hac"] < 0:
        ret_obj.frag_type = "ADDITION"
    else:
        ret_obj.frag_type = "REPLACE"
    return ret_obj


def organise(records, num_picks):
    out_d = {}
    smi_set = set()
    for rec in records:
        rec_key = str(rec)
        addition = {"change": rec.change_frag, "end": rec.end_smi}
        if rec_key in out_d:
            out_d[rec_key]["addition"].append(addition)
        else:
            out_d[rec_key] = {"vector": rec.iso_label, "addition": [addition]}
        smi_set.add(rec.end_smi)
    if num_picks:
        max_per_hypothesis = num_picks / len(out_d)
    out_smi = []
    for rec in out_d:
        # TODO here is the logic as to ordering replacements
        if num_picks:
            random.shuffle(out_d[rec]["addition"])
            out_d[rec]["addition"] = out_d[rec]["addition"][:max_per_hypothesis]
        else:
            out_d[rec]["addition"] = out_d[rec]["addition"]
        out_smi.extend(out_d[rec])
    return out_d


def get_picks(smiles, num_picks, graph_url="neo4j"):
    smiles = canon_input(smiles)
    driver = get_driver(graph_url)
    with driver.session() as session:
        records = []
        for record in session.read_transaction(find_proximal, smiles):
            ans = define_proximal_type(record)
            records.append(ans)
        for record in session.read_transaction(find_double_edge, smiles):
            ans = define_double_edge_type(record)
            records.append(ans)
        for label in list(set([x.label for x in records])):
            # Linkers are meaningless
            if "." in label:
                continue
        if records:
            orga_dict = organise(records, num_picks)
            return orga_dict
        else:
            print("Nothing found for input: " + smiles)


def get_full_graph(smiles, graph_url="neo4j"):
    smiles = canon_input(smiles)
    driver = get_driver(graph_url)
    with driver.session() as session:
        records = []
        for record in session.read_transaction(find_proximal, smiles):
            ans = define_proximal_type(record)
            records.append(ans)
        for record in session.read_transaction(find_double_edge, smiles):
            ans = define_double_edge_type(record)
            records.append(ans)
        for label in list(set([x.label for x in records])):
            # Linkers are meaningless
            if "." in label:
                continue
        if records:
            orga_dict = organise(records, None)
            return orga_dict
        else:
            print("Nothing found for input: " + smiles)


def custom_query(query, graph_url="neo4j"):
    driver = get_driver(graph_url)
    records = []
    with driver.session() as session:
        for record in session.read_transaction(find_custom, query):
            records.append(record)
    return records


def write_picks(smiles, num_picks):
    img_dict = write_results(get_picks(smiles, num_picks))
    for key in img_dict:
        out_f = open(key + ".svg", "w")
        out_f.write(img_dict[key])
