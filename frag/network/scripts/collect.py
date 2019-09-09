import argparse

if __name__ == "__main__":
    # Read in a SD or SMILES file - then write out into a specified directory
    parser = argparse.ArgumentParser(
        description="Annotate a Nodes.txt file with data from multiple attributes files"
    )
    parser.add_argument("--input_attrs")
    parser.add_argument("--input_nodes")
    args = parser.parse_args()
    attrs = []
    id = 0
    # Add multipl annotation readers
    attrs = open(args.input_attrs).readlines()
    attr_dict = {}
    for attr in attrs:
        line_spl = attr.split()
        attr_dict[line_spl[1]] = {"cmpd_id": line_spl[3], "annotation": "EM;MOL;F2"}
    nodes = open(args.input_nodes).readlines()
    out_f = open("nodes.csv", "w")
    out_f.write("smiles:ID(F2),hac:INT,chac:INT,osmiles,cmpd_id,:LABEL\n")
    for x in nodes:
        line_spl = x.split()
        line_new = line_spl[1:]
        if line_new[0] in attr_dict:
            line_new.append(attr_dict[line_new[0]]["cmpd_id"])
            line_new.append(attr_dict[line_new[0]]["annotation"])
        else:
            line_new.append("")
            line_new.append("F2")
        # Now write this out
        out_f.write(",".join(line_new) + "\n")
