from frag.network.decorate import decorate_smi
from frag.utils.vector_utils import get_max_ev_smi

if __name__ == "__main__":
    for item in [
        x.strip().split(",")
        for x in open(
            "/Users/abradley/fragalysis/frag/tests/data/enamine.csv"
        ).readlines()
    ]:
        print(item)
        smiles = decorate_smi(item[1])
        code = item[0]
        for smi in smiles:
            print(smi)
            print((get_max_ev_smi(smi)))
