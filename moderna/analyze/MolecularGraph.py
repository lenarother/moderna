#!/usr/bin/env python
#
# MolecularGraph.py
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Kristian Rother"
__copyright__ = "Copyright 2008, Kristian Rother"
__credits__ = ["Sabrina Hofmann"]
__license__ = "GPL"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"

"""
Topology matching algorithm

(c) 2008 Kristian Rother
GPL Gnu Public License

The topology matcher described here is bronco powerful mojo to deal with 
small molecular structures. 
Among the things you can do are:

- recognize chemical groups in molecules
- recognize modified bases
- annotate chemical attributes in molecules
- read .mol files into a clean data structure


Details:

1. Input data
The mol_parser module reads .mol files with bond information (only the first structure in each file). 
It also can use Residue objects created by the Bio.PDB parser as input. The latter assigns 
single and double bonds using a table of interatomic distances. The resulting data structure 
consists of Molecule (a list of atoms with extras), Atoms (a dictionary with a list 
of bonds), and Bonds (connecting two atoms with a valence).

2. Chemical group recognizer
In the mol_topologies module, there are several hundred topological patterns of chemical 
groups and atom types. They look similar to SMILES, but they can describe bonds with unknown 
valence, and exclude particular groups. It can take care of rings as well.
In the annotated_molecule module, there is a big recursive algorithm that matches the 
topological patterns to atoms within the molecule. The recognition is done in two steps. 
First, the chemical groups are detected in all atoms, and added to the atoms as attributes. 
Second, for all carbons a set of groups in alpha position is examined, resulting in 
additional attributes like 'alpha-carboxy-ketone'.

3. Recognition of modified bases
This is a special set of topological patterns recognizing 120 different RNA base 
modifications. In test_data/molfiles there are sample structures. They originate from 
www.genesilico.pl/modomics, and have been included here as a demonstration of the algorithms 
capabilities.

4. Chemical attributes
This is a set of simple chemical stuff, like electronegativities, oxidative numbers, charges, 
hybridisation etc. They are all added to the atoms as attributes by the annotated_molecule 
module.

5. Unit tests
There are two unit tests that check the most crucial parts of the program. They require the 
PyCogent Python library to work.

6. Usage example:
Getting attributes for a molecule:

from annotated_molecule import AnnotatedMolecule

mol = AnnotatedMolecule()
mol.parse_molfile('test_data/biopath_double.mol')
mol.annotate()
for atom in mol:
    print atom.element, len(atom.bonds)
    print atom['attributes']


"""
import os, re
from MolTopologies import *
from MolParameters import *
from MolGraphParser import *

PERMUTATIONS = (
    [[0]],
    ((0,1),(1,0)),
    ((0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0)),
    ((0,1,2,3),(0,1,3,2),(0,2,1,3),(0,2,3,1),(0,3,1,2),(0,3,2,1),(1,0,2,3),(1,0,3,2),(1,2,0,3),(1,2,3,0),(1,3,0,2),(1,3,2,0),(2,0,1,3),(2,0,3,1),(2,1,0,3),(2,1,3,0),(2,3,0,1),(2,3,1,0),(3,0,1,2),(3,0,2,1),(3,1,0,2),(3,1,2,0),(3,2,0,1),(3,2,1,0)),
    ((0,1,2,3,4),(0,1,2,4,3),(0,1,3,2,4),(0,1,3,4,2),(0,2,1,3,4),(0,2,1,4,3),(0,2,3,1,4),(0,2,3,4,1),(0,3,1,2,4),(0,3,1,4,2),(0,3,4,1,2),(0,3,4,2,1),(0,4,1,2,3),(0,4,1,3,2),(0,4,2,1,3),(0,4,2,3,1),(0,4,3,1,2),(0,4,3,2,1),
     (1,0,4,2,3),(1,0,4,3,2),(1,0,2,4,3),(1,0,2,3,4),(1,0,3,4,2),(1,0,3,2,4),(1,4,0,2,3),(1,4,0,3,2),(1,4,2,0,3),(1,4,2,3,0),(1,4,3,0,2),(1,4,3,2,0),(1,2,0,4,3),(1,2,0,3,4),(1,2,4,0,3),(1,2,4,3,0),(1,2,3,0,4),(1,2,3,4,0),(1,3,0,4,2),(1,3,0,2,4),(1,3,4,0,2),(1,3,4,2,0),(1,3,2,0,4),(1,3,2,4,0),
     (2,0,4,1,3),(2,0,4,3,1),(2,0,1,4,3),(2,0,1,3,4),(2,0,3,4,1),(2,0,3,1,4),(2,4,0,1,3),(2,4,0,3,1),(2,4,1,0,3),(2,4,1,3,0),(2,4,3,0,1),(2,4,3,1,0),(2,1,0,4,3),(2,1,0,3,4),(2,1,4,0,3),(2,1,4,3,0),(2,1,3,0,4),(2,1,3,4,0),(2,3,0,4,1),(2,3,0,1,4),(2,3,4,0,1),(2,3,4,1,0),(2,3,1,0,4),(2,3,1,4,0),
     (3,0,4,2,1),(3,0,4,1,2),(3,0,2,4,1),(3,0,2,1,4),(3,0,1,4,2),(3,0,1,2,4),(3,4,0,2,1),(3,4,0,1,2),(3,4,2,0,1),(3,4,2,1,0),(3,4,1,0,2),(3,4,1,2,0),(3,2,0,4,1),(3,2,0,1,4),(3,2,4,0,1),(3,2,4,1,0),(3,2,1,0,4),(3,2,1,4,0),(3,1,0,4,2),(3,1,0,2,4),(3,1,4,0,2),(3,1,4,2,0),(3,1,2,0,4),(3,1,2,4,0),
     (4,0,1,2,3),(4,0,1,3,2),(4,0,2,1,3),(4,0,2,3,1),(4,0,3,1,2),(4,0,3,2,1),(4,1,0,2,3),(4,1,0,3,2),(4,1,2,0,3),(4,1,2,3,0),(4,1,3,0,2),(4,1,3,2,0),(4,2,0,1,3),(4,2,0,3,1),(4,2,1,0,3),(4,2,1,3,0),(4,2,3,0,1),(4,2,3,1,0),(4,3,0,1,2),(4,3,0,2,1),(4,3,1,0,2),(4,3,1,2,0),(4,3,2,0,1),(4,3,2,1,0))
    )

    
class AnnotatedAtom(Atom):

    def nb_has_feature(self,feature):
        """Auxiliary procedure for condition tree."""
        for bond in self.bonds:
            if feature in bond.atom2['attributes']:
                return True

    def calc_charge(self):
        """
        Count electrons that participate in bonds, and compare to the number
        of outer electron each atom should have, add the 
        'charge' and 'free_electron_pairs' fields.
        Radicals are not considered here
        """
        if OUTER_ELECTRONS.has_key(self.element):
            outer = OUTER_ELECTRONS[self.element]
        else:
            outer = 1
            print "UNKNOWN ELEMENT",self.element
        electrons  = 0
        free_pairs = 0
        charge     = 0
        for bond in self.bonds:
            electrons += 2*bond.valence

        if self.element == 'H':
            if electrons == 1: charge = +1
            # elif electrons != 2: print 'WARNING: WEIRD HYDROGEN CONFIGURATION!'
        elif self.element == 'C':
            if electrons == 6: charge = +1
        elif self.element == 'N':
            if electrons == 6: free_pairs = 1
            elif electrons == 8: charge = +1
        elif self.element == 'O':
            if electrons == 4: free_pairs = 2
            elif electrons == 2:
                free_pairs = 3
                charge     = -1
        elif self.element == 'S':
            if electrons == 4: free_pairs = 2            
        elif self.element in ['Fe','Mg','Ca','K','Li','Na']: 
            charge = +2

        # charges &  free electron pairs
        self['charge'] = charge
        self['free_electron_pairs'] = free_pairs

    def calc_electronegativity(self):
        """assign the electronegativity field."""
        self['electronegativity'] =  ELECTRONEGATIVITIES.get(self.element, 1.0)

    def calc_oxidative_number(self):
        """
        Sum up the electronegativity differences for the neighbors of each bond
        to calculate the 'partial charge' and 'oxidative number' attributes.
        """
        self['partial_charge'] = 0.0
        # oxidative number
        oxi = self['charge']
        for bond in self.bonds:
            # compare electronegativities
            bonded_en = bond.atom2['electronegativity']
            en_diff = bonded_en - self['electronegativity']
            if en_diff < 0 : 
                oxi -= bond.valence
            elif en_diff > 0 : 
                oxi += bond.valence
            # locate polar bonds.
            # according to wikipedia, polar bonds have en-differences
            # between 0.5 and 1.7.
            # For ion bonds (en-diff>1.7), a warning will be printed
            if abs(en_diff) >= POLAR_BOND_MIN:
                if abs(en_diff) > POLAR_BOND_MAX: 
                    print "WARNING: Ion bond marked as bond with e.n.-difference of",en_diff
                else: 
                    self['partial_charge'] += en_diff*bond.valence
                    
        self['oxidative_number'] = oxi
        # round partial charge to 1 digit
        self['partial_charge'] = round(self['partial_charge'],1)

    def calc_hybridisation(self):
        """
        Assigns the 'hybridisation' attribute.
        """
        # count the bonds of each kind
        bonds = [0,0,0,0,0,0]
        for bond in self.bonds:
            if bond.valence<len(bonds):
                bonds[bond.valence] += 1
            else:
                print "HUGE VALENCE FOUND",bond.valence
        total_bonds = bonds[1] + bonds[2]*2 + bonds[3]*3
        
        # check for single atoms
        if total_bonds == 0 and self.element not in ['C','O','N']:
            self.add_attribute('single atom')
            if self['charge']!=0: 
                self.add_attribute('ion')
        # assign hybridisation state
        if self.element not in ['H','He']:
            if total_bonds > 4: 
                self.add_attribute('higher_order_hybrid')
            elif bonds[3] > 0 or bonds[2]==2: 
                self.add_attribute('sp-hybrid')
            elif bonds[2] == 1:
                self.add_attribute('sp2-hybrid')
            else: 
                self.add_attribute('sp3-hybrid')
            
    def detect_h_don_acc(self):
        """
        Assigns H-donor and  H-acceptor attributes.
        """
        # identify h donors
        if ELECTRONEGATIVITIES.has_key(self.element):
            diff = ELECTRONEGATIVITIES[self.element]-ELECTRONEGATIVITIES['H']
        else:
            diff = -1
        
        if diff >= H_DONOR_DIFFERENCE:
            donor = 0
            for bond in self.bonds:
                if bond.atom2.element == 'H':
                    donor = 1
                    bond.atom2.add_attribute('h-bond hydrogen')
            if donor: 
                self.add_attribute('h-bond donor')
    
        # identify h acceptors
        if diff >= H_ACCEPTOR_DIFFERENCE:
            if self['free_electron_pairs']>0:
                self.add_attribute('h-bond acceptor')


    def detect_chirality(self):
        """
        Finds out if the atom has four different substituents.
        """
        if len(self.bonds) >= 4:
            # potential center of chirality
            branches = {}
            for bond in self.bonds:
                molstring = bond.atom2.get_molstring([self])
                branches[molstring] = True
            if len(branches.keys())>= 4:
                # four different branches.. chiral.
                self.add_attribute('chiral_atom')

    def is_maybe_conjugated(self):
        """
        Checks whether atom j has a double bond, charge or free electron pair.
        """
        if self['free_electron_pairs'] > 0: return True
        if self['charge'] != 0: return True
        for bond in self.bonds:
            if bond.valence >= 2: return True
        return False

    def get_conjugated_system(self, candidates,visited,depth = -1):
        """
        Returns a list of atoms belonging to a conjugated system
        """
        if depth == 0: 
            return []
        if self in visited: 
            return []
        conj_sys = [self]
        for bond in self.bonds:
            if bond.atom2 in candidates and bond.atom2 not in visited:
                more_atoms = bond.atom2.get_conjugated_system(candidates,visited[:]+[self],depth-1)
                for m in more_atoms:
                    if m not in conj_sys: 
                        conj_sys.append(m)
        return conj_sys

    def match_node(self,node,visited,indices):
        """
        Graph matching algorithm for detecting functional groups.
        The functional groups are given as starting nodes in a
        topology of bonded atoms, parsed into a data structure
        by parameters.py.
        
        This recursive procedure checks the neighbors at each node
        one by one. returns 0 or 1.

        CAVEAT: This algorithm cannot distinguish whether two
        sub-branches of a structure meet not. E.g. the topologies
        C(-C(-C3),-C(-C3)) and C(-C(-C),-C(-C)) cannot be distinguished.
        Instead, you should write the first one as C3(-C(-C(-C(-C3)))),
        the second one as C(-C(-C(-C(-C)))).
        """
        # The node is a parsed version of the patterns in the mol_topologies module
        # Each node contains [element, list of subnodes, negation bit, index]
        elem, subnodes, boolean, index, valence = node
        # 1. check the element
        if elem != self.element: return 0
        
        # 2. check index number
        # the 'indices' field contains candidates for atoms that have
        # a number in the pattern, e.g. C1(-C(-C(-C(-C(-C1))))) for cyclopentan.
        if index > 0:
            for idx_index, idx_atom  in indices:
                ia = idx_atom == self
                ib = idx_index == index
                if ia or ib:
                    if ia and ib: return 1
                    return 0
            if self in visited:
                return 0
            else:
                # this is a new indexed atom
                indices.append((index,self))            
        
        # 3. check the subnodes of this atom
        # it needs to be ensured that for each subnode,
        # a different bond is used. Otherwise, a single
        # methyl group would be enough to fulfill a (-C,-C,-C) pattern

        # shortcut for optimized performance
        if len(subnodes) == 0: return 1

        # now, each of the nodes needs to match a different bond.
        # This is done by building a [nodes][bonds] binary matrix
        # of which nodes match which bonds.
        # Additionally, a list of the excluded bonds is created.
        bmatrix = []
        excluded_bonds = [0]*self.n_bonds
        include_nodes = []

        for subnode in subnodes:
            sub_elem, sub_subnodes, sub_boolean, sub_index,sub_valence = subnode
            if sub_boolean:
                # this is an included node that must be matched
                include_nodes.append(subnode)
                row = []
                for bond in self.bonds:
                    # first condition is valence.if it doesnt match no further checks are needed
                    if sub_valence<=0 or bond.valence==sub_valence:
                        # if this atom was already visited, it must be indexed
                        if bond.atom2 in visited:
                            if sub_index:
                                # check the indexed atoms
                                row.append(bond.atom2.match_node(subnode,visited[:]+[self],indices[:]))
                            else:
                                # unindexed visited bonds never match!
                                row.append(0)
                        else:
                            # check unvisited bonds, no matter whether they are indexed or not
                            row.append(bond.atom2.match_node(subnode,visited[:]+[self],indices[:]))
                    else:
                        # bonds with a wrong valence dont match
                        row.append(0)

                # shortcut: if no match at all in this node, forget it.
                if sum(row) == 0: return 0
                bmatrix.append(row)
                
            else:
                # this is an excluded node that must not be matched.
                # update list of bonds that match excluded nodes
                j = 0
                for bond in self.bonds:
                    # first condition is valence, again
                    if sub_valence<=0 or bond.valence==sub_valence:
                        # if this atom was already visited, it could be indexed
                        if bond.atom2 in visited:
                            if sub_index:
                                # check the indexed atoms. if an excluded indexed node, e.g. !-C3 is  
                                # found again, this is marked in the exclusion list. = BAD
                                if bond.atom2.match_node(subnode,visited[:]+[self],indices[:]):
                                    excluded_bonds[j] = 1
                        elif bond.atom2.match_node(subnode,visited[:]+[self],indices[:]):
                            # unvisited bond matches something excluded = BAD
                            excluded_bonds[j] = 1
                    j += 1

        # now check the binary matrix. The condition is, that for every node 
        # a True must exist in a different column. Dont know what this is called in math lingo.
        # Practically, each permutation of the columns is looped,
        # and a True diagonal from [0][0] to [n][n] is looked for.
        if len(include_nodes) > self.n_bonds: return 0
        if len(self.bonds)>5:
            print "MORE THAN FIVE BONDS ON ATOM, SKIPPING."
            return 0
        for permu in PERMUTATIONS[self.n_bonds-1]:
            diagonal = 1
            for j in range(len(include_nodes)):
                if not bmatrix[j][permu[j]]: diagonal = 0
            if diagonal:
                # now check whether there any of the remaining bonds
                # have been excluded
                for j in range(len(include_nodes),self.n_bonds):
                    if excluded_bonds[permu[j]]:
                        # there was one bad group remaining that spoils everything
                        return 0
                return 1
        return 0

        
class AnnotatedMolecule(Molecule):
    atom_class = AnnotatedAtom
    
    def __init__(self):
        Molecule.__init__(self)
        self.rings = []
        self.sum_formula = {}
    
    def detect_conjugated(self):
        """conjugated bonds/mesomeric systems are assigned, where an atom
        - has a double bond, charge or free electron pair
        AND
        - the same condition applies to one of its single-bonded neighbors
        """
        # detect conjugated
        system_number = 1
        candidates = [] 
        for atom in self:
            if atom.is_maybe_conjugated(): 
                candidates.append(atom)

        self.conjugated_systems = []
        # now group the candidates to connected systems with 3+ atoms
        while len(candidates)>2:
            csystem = candidates[0].get_conjugated_system(candidates,[])
            if len(csystem)>=3:
                self.conjugated_systems.append(csystem)
                for atom in csystem:
                    if 'conjugated_system' not in atom['attributes']:
                        atom.add_attribute('conjugated_system')

            # now pick remaining candidates
            new_cand = []
            for c in candidates:
                if c not in csystem: 
                    new_cand.append(c)
            candidates = new_cand
            
    def get_sum_formula(self):
        for a in self:
            self.sum_formula.setdefault(a.element,0)
            self.sum_formula[a.element] += 1

    def check_sum_formula(self, sumdict):
        """Makes sure that at least the same number of atoms
        as in the topologies sum formula are present.
        Returns boolean
        """
        for s in sumdict:
            if s == 'H': continue
            if self.sum_formula.get(s,0) < sumdict[s]: return False
        for s in ['S','Se']:
            if self.sum_formula.has_key(s) and not sumdict.has_key(s):
                return False
        return True
        
    def detect_patterns(self, pattern_list, c_only=False):
        """
        Detects matching patterns in a list of Pattern objects.
        Calls the topology matching algorithm on an atom for each topology  in topologies.
        pattern_list - list of parsed topologies from mol_topologies.
        c_only - only consider topologies starting at carbons.
        """
        result = []
        self.get_sum_formula()
        for pattern in pattern_list:
            # for accelerated search: check sum formula first
            if not self.check_sum_formula(pattern.sumdict): continue
            for i, atom in enumerate(self):
                if not c_only or atom.element == 'C':
                    if atom.match_node(pattern.topology,[],[]):
                        #atom.add_attribute(pattern.name)
                        result.append((pattern.name, i))
        return result

    def detect_functional_groups(self):
        """
        Runs the topology matching algorithm on functional groups, by combining 
        each functional group topology with each topology for groups in alpha position.
        """
        functional = self.detect_patterns(FUNCTIONAL_GROUP_TOPOLOGIES)
        alpha = self.detect_patterns(ALPHA_TOPOLOGIES,True)
        # add normal group attributes
        for func, i in functional:
            atom = self[i]
            atom.add_attribute(func)
            # go through atoms in alpha position for all carbons
            if atom.element == 'C':
                for bond in atom.bonds:
                    j_bond = self.index(bond.atom2)
                    for alphagroup, j in alpha:
                        if j == j_bond:
                            name = 'alpha-'+alphagroup+'-'+func
                            atom.add_attribute(name)

    def traverse_rings(self,rings,traceback,dead_ends,depth=-1):
        """
        Detects rings in the molecules
        Recursively swoops along the bonds for loops.
        Anytime the algorithm bites its tail a ring is found.
        """
        if depth == 0: return
        atom = traceback[-1]
    
        if len(atom.bonds)<=1:
            # any atom with one or no neighbor cant be in a ring
            dead_ends.append(traceback[-1])
        else:            
            # travel all bonds    
            for bond in atom.bonds:
                if len(traceback)>2 and bond.atom2 == traceback[0]:
                    # whopee! a ring has been found!
                    ring = traceback[:]
                    ring.sort() # sortiert die atom-indices
                    string_ring = string.join([str(rr) for rr in ring],',')
                    rings[string_ring] = traceback[:]
                elif bond.atom2 in dead_ends: continue
                elif bond.atom2 in traceback: continue                
                else:
                    # travel on
                    self.traverse_rings(rings,traceback[:]+[bond.atom2],dead_ends,depth-1)

            
    def detect_rings(self):
        # detect ring structures
        # each atom is theoretically a starting point for rings.
        dead_ends = []
        rings = {}
        for atom in self:
            self.traverse_rings(rings,[atom],dead_ends,MAX_RING_SIZE)

        # match rings to patterns
        # a fingerprint for each ring is created,
        # e.g. C=C-C=C-C=C- for a benzene ring.
        for r in rings.keys():
            ring = rings[r]
            ring_size = len(ring)
            ring_fingerprint = ""

            # traverse the ring and check the bonds
            for i in range(ring_size):
                ring_atom = ring[i]
                ring_fingerprint += ring_atom.element
                next = i+1
                if next>=ring_size: next=0
                for bond in ring_atom.bonds:
                    if bond.atom2 == ring[next]:                        
                        if bond.valence==1: bstr = '-'
                        elif bond.valence==2: bstr ='='
                        else: bstr='#'
                        ring_fingerprint += bstr
            # if any of the according bonds are not found,
            # the detection algorithm must be buggy
            if len(ring_fingerprint) != 2*ring_size:
                print "ERROR IN RING DETECTION"
                print ring_fingerprint
                print ring
            
            # do a circular permutation of the fingerprint
            # to match it with the patterns
            i = 0
            matched = 0
            while i<ring_size and not matched:
                if RING_PATTERNS.has_key(ring_fingerprint):
                    matched = RING_PATTERNS[ring_fingerprint]
                else:
                    ring_fingerprint = ring_fingerprint[-2:]+ring_fingerprint[:-2]
                    i += 1

            # store the ring
            if not matched: matched = 'unknown_ring'
            self.rings.append((ring,ring_fingerprint,matched))
        
            # annotate atoms
            aromatic = 1 # rings are aromatic unless concluded otherwise
            for atom in ring:
                if not 'part_of_ring' in atom['attributes']:
                    atom['attributes'].append('part_of_ring')
                if matched:
                    atom['attributes'].append(matched)

                # check whether it is an aromatic ring
                if not 'conjugated_system' in atom['attributes']:
                    aromatic = 0

            if aromatic:
                for atom in ring:
                    if 'aromatic_ring' not in atom['attributes']:
                        atom['attributes'].append("aromatic_ring")


    def test_condition(self,condition,atom,elem,bondtype,attributes,neighbors):
        result = 0
        exec "if %s: result = 1"%(condition)
        return result

    
    def assign_wang_atomtype(self):
        # hierarchical atomtyping schema.
        #
        # The assignment is according to the table in Wang et al. 2004
        # 'Development and testing of a general Amber force field'
        # 
        # It is done by using Python to generate and execute python code
        # from the conditions listed above.
        #
        for atom in self:
            elem = atom.element
            bondtype = atom['bondschema']
            attributes = atom['attributes']
            neighbors = []  
            for bond in atom.bonds:
                if bond.valence==1: bond_order='-'
                elif bond.valence==2: bond_order='='
                elif bond.valence==3: bond_order='#'
                else: bond_order='~'
                neighbors.append(bond_order+bond.atom2.element)
            
            def recurse_tree(tree):
                """Recursive procedure to process the tree above and to assign atom types."""
                first = tree[0]
                last = tree[-1]
                ic = 0
                if type(tree)==type("abc"):
                    return tree
                elif type(first)==type("abc"):
                    # test a condition from the tree
                    if self.test_condition(first,atom,elem,bondtype,attributes,neighbors):
                        att = None
                        return recurse_tree(last)
                    else: return None
                else:
                    att = None
                    for condition in tree[:-1]:
                        if self.test_condition(condition[0],atom,elem,bondtype,attributes,neighbors):
                            att = recurse_tree(condition[1])
                            break
                    if not att: att = last
                    return att
            
            # print elem,bondtype,neighbors,attributes
            atomtype = recurse_tree(TYPE_DECISION_TREE)
            atom['wang_atomtype'] = atomtype
                
                
    def annotate(self):
        for atom in self: atom.calc_charge()
        for atom in self: atom.calc_electronegativity()
        for atom in self: atom.calc_oxidative_number()    
        for atom in self: atom.calc_hybridisation()        
        for atom in self: atom.detect_bonding_schema()     
        for atom in self: atom.detect_h_don_acc()          
        for atom in self: atom.detect_chirality()
        self.detect_conjugated()         
        self.detect_functional_groups()  
        self.detect_rings()              
        self.assign_wang_atomtype()


#
#
# EASTER EGG SECTION
#
#
def is_subring_of(ring1,ring2):
    """Checks whether ring1 is fully contained in ring2)"""
    subring = 1
    for r1 in ring1:
        atom = r1[0]
        found = 0
        for r2 in ring2:
            if r2[0]==atom: found = 1
        if not found: subring = 0
    return subring

def complexrings():
    #
    # rogue code
    # 
    # detect complicated ring structures
    for r2 in rings.keys():
        subrings = []
        # look for sub-rings
        for r1 in rings.keys():
            if is_subring_of(r1,r2): subrings.append(rings[r1][2])
        if rings[r2][2]=='putative_purine':
            if "pyrimidine" in subrings and "5-ring of purine" in subrings:
                rings[r2][2]='purine'
        elif rings[r2][2]=='putative_flavin':
            if "benzene" in subrings and "putative_central_flavin_ring" in subrings and "putative_two_rings_of_flavin" in subrings and "putative_pyrimidine_ring_of_flavin" in subrings:
                rings[r2][2]='isoalloxazin'
        elif rings[r2][2]=='putative_tetrapyrrol':
            pyrrols = []
            for r1 in rings.keys():
                if rings[r1][2]=='pyrrol':
                    for i in rings[r1][0]: pyrrols.append(i[0])
            pyrrol = 0
            aromatic = 1
            for at in rings[r2][0]:
                if at[0] in pyrrols: pyrrol += 1
                if not has_feature(mol[at[0]],"aromatic_ring"): aromatic = 0
            if pyrrol == 16 and aromatic: rings[r2][2] = "tetrapyrrol"
            elif pyrrol == 16:
                print rings[r2]
        

def find_patterns(mol):
    """
    Detects and annotates all kinds of chemical group.
    The procedure goes several times over the entire molecule,
    first annotating the easy features, then the difficult ones.
    """
    # classify some easy groups
    for atom in mol:
        typ = atom[7]
        if ASSOCIATIONS.has_key(typ):
            add_feature(atom,ASSOCIATIONS[typ])
                             
    # do some statistics
    for atom in mol:
        for p in string.split(atom[6],','):
            if not all_features.has_key(p): all_features[p] = 0
            all_features[p] += 1


if __name__ == '__main__':
    print """
    Annotating chemical groups in mol files
    usage: 
    annotate_molecule.py <.mol or .mol2 file> 
    or 
    annotate_molecule.py <path with many mol files>
    """
    if len(sys.argv)>1:
        mol = AnnotatedMolecule()
        if sys.argv[1][-4:] == '.mol':
            mol.parse_molfile(sys.argv[1])
        elif sys.argv[1][-5:] == '.mol2':
            mol.parse_mol2(sys.argv[1])
        mol.annotate()
        for atom in mol:
            print atom, atom['attributes']
    else:
        for f in os.listdir(sys.argv[1]):
            print f
            mol = Molecule()
            if f[-4:] == '.mol':
                mol.parse_molfile('mol'+os.sep+f)
            elif f[-5:] == '.mol2':
                mol.parse_mol2('mol'+os.sep+f)
            mol.annotate()
            mol.report()

        for k in all_types.keys():print k,all_types[k]
        for f in all_features.keys():print f,all_features[f]

