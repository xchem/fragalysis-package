from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def get_best_fit_plane(pts, weights=None):
    """

    :param pts:
    :param weights:
    :return:
    """
    linalg = np.linalg
    if weights is None:
        wSum = len(pts)
        origin = np.sum(pts, 0)
    origin /= wSum
    sums = np.zeros((3, 3), np.double)
    for pt in pts:
        dp = pt - origin
        for i in range(3):
            sums[i, i] += dp[i] * dp[i]
            for j in range(i + 1, 3):
                sums[i, j] += dp[i] * dp[j]
                sums[j, i] += dp[i] * dp[j]
    sums /= wSum
    vals, vects = linalg.eigh(sums)
    order = np.argsort(vals)
    normal = vects[:, order[0]]
    plane = np.zeros((4,), np.double)
    plane[:3] = normal
    plane[3] = -1 * normal.dot(origin)
    return plane


def get_pbf(plane, pts):
    denom = np.dot(plane[:3], plane[:3])
    denom = denom ** 0.5
    # add up the distance from the plane for each point:
    res = 0.0
    for pt in pts:
        res += np.abs(pt.dot(plane[:3]) + plane[3])
        res /= denom
        res /= len(pts)
    return res


def plane_best_fit_exit_vector(scaffold, repl_smarts, confId):
    """
    Get the exit vector for the best fit plane and the orthogonal plane.
    :param scaffold:
    :param repl_smarts:
    :param confId:
    :return:
    """
    # Get PBF plane for murcko scaffold only
    conf = scaffold.GetConformer(confId)
    if not conf.Is3D():
        print("This mol is not 3D - all PBFev angles will be 0 degrees")
        return [0]
    pts = np.array(
        [
            list(conf.GetAtomPosition(i))  # Get atom coordinates if not Xenon
            for i in range(scaffold.GetNumAtoms())
            if scaffold.GetAtomWithIdx(i).GetSymbol() != repl_smarts
        ]
    )
    # Plane is xyz vector with a c intercept adjustment
    plane = get_best_fit_plane(pts)
    pbf = get_pbf(plane, pts)
    # Where [#0] matches exit vector SMILES [*]
    patt = Chem.MolFromSmarts("[" + repl_smarts + "]-[*]")
    matches = scaffold.GetSubstructMatches(patt)
    if len(matches) == 0:
        return None

    # Calculate angles between exit vectors and the murcko plane of best fit
    exitVectors = {}
    denom = np.dot(plane[:3], plane[:3])
    denom = denom ** 0.5
    for n, match in enumerate(matches):
        out_isotope = scaffold.GetAtomWithIdx(match[0]).GetIsotope()
        evCoords = conf.GetAtomPosition(match[0])
        anchorCoords = conf.GetAtomPosition(match[1])
        v = np.array(
            (
                (evCoords[0] - anchorCoords[0]),
                (evCoords[1] - anchorCoords[1]),
                (evCoords[2] - anchorCoords[2]),
            )
        )
        angle = np.arcsin((np.dot(v, plane[:3])) / ((denom) * ((np.dot(v, v)) ** 0.5)))
        angle = np.abs(np.degrees(angle))
        exitVectors[out_isotope] = angle
    exitVectors["pbf"] = pbf
    return exitVectors


def get_exit_vector_for_xe_smi(input_smi, repl_smarts="At"):
    """
    Get the exit vector angles for a SMILES with Xen
    :param input_smi:
    :return:
    """
    mol = Chem.MolFromSmiles(input_smi.replace("Xe", repl_smarts))
    conf_id = AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
    return plane_best_fit_exit_vector(mol, repl_smarts, conf_id)


def get_max_ev_smi(smi):
    """
    Try three different size groups - and take the maximum.
    This is a fairly huge assumption. The variability in the vector coud also be a parameter
    :param smi:
    :return:
    """
    # Big group
    big = get_exit_vector_for_xe_smi(smi.replace("At", "Xe"), "At")
    # Mediume group
    medium = get_exit_vector_for_xe_smi(smi.replace("At", "Xe"), "K")
    # Little group
    small = get_exit_vector_for_xe_smi(smi.replace("At", "Xe"), "Li")
    out_d = {}
    for key in big:
        out_d[key] = max([big[key], medium[key], small[key]])
    return out_d
