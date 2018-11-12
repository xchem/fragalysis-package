import json

from frag.utils.network_utils import write_results


def write_all_out(out_d):
    for key in out_d:
        img = write_results(out_d[key])
        out_file = open(key+".svg","w")
        out_file.write(img)

if __name__ == "__main__":
    write_all_out(json.load(open("data.json")))

