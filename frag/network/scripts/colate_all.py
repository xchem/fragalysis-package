import os
import shutil

def do_for_dir(input_dir):
    os.chdir(input_dir)
    prop_dict = {"nodes.txt":[None,"smiles:ID(F2)", "hac:INT", "chac:INT", "osmiles"],
    "edges.txt":[None,":START_ID(F2)",":END_ID(F2)","label"],
    }
    attr_smis = []
    attr_ids = []
    for x in open("attributes.txt").readlines():
        attr_smis.append(x.split()[1])
        attr_ids.append(x.split()[3])

    for f_name in prop_dict:
        out_f = open(f_name.replace(".txt",".csv"),"w")
        for line in open(f_name).readlines():
            line_spl = line.split()
            out_l = []
            for i,x in enumerate(prop_dict[f_name]):
                if x:
                    out_l.append(line_spl[i])
            if f_name == "nodes.txt":
                if line_spl[1] in attr_smis:
                    cmpd_id = attr_ids[attr_smis.index(line_spl[1])]
                    out_l.append(cmpd_id)
                    # This is where we can add tags - like CHEAP - EXPENSIVE
                    out_l.append("EM;MOL;F2")
                else:
                    out_l.append("")
                    out_l.append("F2")
                out_f.write(",".join(out_l)+"\n")
            elif f_name =="edges.txt":
                out_f.write(",".join(out_l)+"\n")
        out_f.flush()
        out_f.close()

node_list = ["nodes-header.csv"]
edge_list = ["edges-header.csv"]
with open("edges-header.csv","w") as out_f:
    out_f.write(",".join([x for x in [":START_ID(F2)",":END_ID(F2)","label"]]))
with open("nodes-header.csv", "w") as out_f:
    out_f.write(",".join([x for x in ["smiles:ID(F2)", "hac:INT", "chac:INT",
                                      "osmiles", "cmpd_id", ":LABEL"] if x]) + "\n")
do_for_dir(".")
node_list.append("nodes.csv")
edge_list.append("edges.csv")
# Now this is the command to run
print(" ".join(["/var/lib/neo4j/bin/neo4j-admin import","--database", "new.db", "--nodes",'"'+",".join(node_list)+'"',
                "--relationships:F2EDGE",'"'+",".join(edge_list)+'"']))
