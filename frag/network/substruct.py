import pickle, time
from rdkit.Chem import rdSubstructLibrary
from rdkit import Chem
import tqdm

def generate_fps(input_smi):
    # start by building the fingerprint cache
    t1 = time.time()
    with open(input_smi, 'r') as inf:
        ls = [[x.split()[0],x.split()[1]] for x in inf]
        ls.pop(0)
        with open(input_smi.replace(".smi",".pkl"), 'wb+') as pklf:
            for i, (smi,chemblid) in enumerate(ls):
                m = Chem.MolFromSmiles(smi, sanitize=False)
                fp = Chem.PatternFingerprint(m)
                pickle.dump(fp, pklf)
                if not (i + 1) % 10000:
                    print("Done", i + 1)
    t2=time.time()
    print("That took %.2f seconds. The output has %d molecules."%(t2-t1,i))



def read_in_lib(input_smi):
    t1=time.time()
    mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps = rdSubstructLibrary.PatternHolder()
    with open(input_smi,'r') as inf:
        ls = [x.split() for x in inf]
        ls.pop(0)
        with open(input_smi.replace(".smi",".pkl"),'rb') as pklf:
            for l in tqdm.tqdm(ls):
                smi = l[1]
                mols.AddSmiles(smi)
                fp = pickle.load(pklf)
                fps.AddFingerprint(fp)
    library = rdSubstructLibrary.SubstructLibrary(mols,fps)
    t2=time.time()
    print("That took %.2f seconds. The library has %d molecules."%(t2-t1,len(library)))
    return library
