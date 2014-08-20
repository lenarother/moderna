
# analyzes whether the backbones of RNA residues are connected.

from moderna.RNAResidue import BACKBONE_DIST_MATRIX, DIST_TOLERANCE, O3_P_DIST_HI

def are_residues_connected(res1, res2):
    """
    Checks whether two residues are connected.
    Distances on the backbone are within norm + tolerance.

    Arguments:
    * upstream residue as RNAResidue object
    * downstream residue as RNAResidue object
    """
    try:
        if res1["O3'"] - res2["P"] > O3_P_DIST_HI * DIST_TOLERANCE:
            return False
        for tup in [("C3'", "O3'"), ("C4'", "C3'")]:
            a1, a2 = tup
            low, hi = BACKBONE_DIST_MATRIX[tup]
            if res1[a1] - res1[a2] > hi * DIST_TOLERANCE:
                return False
        for tup in [   ("P", "O5'"), ("O5'", "C5'"), ("C5'", "C4'")]:
            a1, a2 = tup
            low, hi = BACKBONE_DIST_MATRIX[tup]
            if res2[a1] - res2[a2] > hi * DIST_TOLERANCE:
                return False
        return True
    except KeyError:
        # missing atoms
        return False


def is_chain_continuous(chain):
    """
    Checks whether a chain is continuous.
    Check whether all subsequent pairs of residues are connected.

    Arguments:
        chain - RNAChain object
    """
    keys_sorted = []
    for resi in chain:
        keys_sorted.append(resi.identifier)

    pairs_to_check = [(keys_sorted[x], keys_sorted[x+1]) for x in range(len(keys_sorted)-1)]

    for pair in pairs_to_check:
        if not are_residues_connected(chain[pair[0]], chain[pair[1]]):
            return False
    return True
