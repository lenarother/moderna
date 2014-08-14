#
# topology_graph.py
#
# helper functions for the ModernaGUI
#
from mol_parser import *

class BondTopology(list):
    """represents a bond between two Bio.PDB.Atoms."""
    def __init__(self, a1, a2):
        self.atom1 = a1
        self.atom2 = a2
        self.append(a1.coord)
        self.append(a2.coord)
        

def get_resi_coordinates(struc):
    """Analyzes bonds in each residue."""
    result = []
    for resi in struc:
        m = Molecule()
        m.parse_resi(resi)
        # get all bonds
        all_bonds = []
        for a in m:
            for b in a.bonds:
                # remove duplicate bonds
                duplicate = False
                if b in all_bonds: duplicate = True
                else:
                    for bond in all_bonds:
                        if b.atom1 == bond.atom2 and b.atom2==bond.atom1:
                            duplicate = True
                if not duplicate:
                    all_bonds.append(b)
        # associate atoms with atoms from bond recognizer
        for pdb_at,  mol_at in zip(resi, m):
            # replace atoms in bonds by pdb atoms
            # THIS IS UGLY OOC!!
            for b in all_bonds:
                if b.atom1 == mol_at:
                    b.atom1 = pdb_at
                elif b.atom2 == mol_at:
                    b.atom2 = pdb_at
        # create result
        for b in all_bonds:
            result.append(BondTopology(b.atom1, b.atom2))
    return result 

def get_interresidue_bonds(struc):
    """Analyzes bonds between two residues."""
    result = []
    for resi1 in struc:
        for resi2 in struc[resi1.identifier:]:
            if resi1==resi2: continue
            if struc.are_residues_connected(resi1, resi2):
                try:
                    result.append(BondTopology(resi1["O3'"], resi2["P"]))
                except KeyError:
                    # no aproppriate atoms in residues --> no bond
                    pass
    return result
        
        
def get_bond_coordinates(struc):
    """Creates a list of BondToplogy objects from a ModernaStructure."""
    result = []
    result += get_resi_coordinates(struc)
    result += get_interresidue_bonds(struc)
    return result
    
