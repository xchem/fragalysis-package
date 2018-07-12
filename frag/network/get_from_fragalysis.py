import requests

prot_code = "MURD-x0373"
r = requests.get(
    "https://fragalysis.apps.xchem.diamond.ac.uk/api/proteins/",
    params={"code": prot_code},
)
result = r.json()["results"][0]
pdb_url = result["pdb_info"]
pdb_id = result["id"]
r = requests.get(
    "https://fragalysis.apps.xchem.diamond.ac.uk/api/molecules/",
    params={"prot_id": pdb_id},
)
out_json = r.json()["results"][0]
sdf_info = out_json["sdf_info"]
out_f = open(prot_code + ".mol", "w")
out_f.write(sdf_info)
out_f.close()
# Now s
r = requests.get(pdb_url)
out_f = open(prot_code + "_apo.pdb", "w")
out_f.write(r.text)
out_f.close()
