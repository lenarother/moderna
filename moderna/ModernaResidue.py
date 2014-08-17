#!/usr/bin/env python
#
# ModernaResidue.py
#
# Residue object with modeling functionality.
#
# http://iimcb.genesilico.pl/moderna/ 
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

from RNAResidue import RNAResidue
from Bio.PDB import PDBParser
from ModernaSuperimposer import ModernaSuperimposer
from sequence.ModernaAlphabet import alphabet
from Bio.PDB.Vector import rotaxis2m
from Errors import ModernaResidueError
import math

from Constants import UNKNOWN_RESIDUE_SHORT, \
    ADDING_MODIFICATION_RULES_PATH, MODIFICATION_FRAGMENTS_PATH, \
    BACKBONE_RIBOSE_ATOMS, ANY_RESIDUE, MISSING_RESIDUE, RIBOSE, \
    BASE_PATH, B_FACTOR_REMOVE_MODIF, B_FACTOR_ADD_MODIF, B_FACTOR_EXCHANGE


class ModernaResidue(RNAResidue):
    """
Deals with residue objects.
Supplements RNAResidue object with functions to
- exchange bases
- add modifications
- remove modifications
- rotat chi angle
The type of each ModernaResidue is automatically recognized, and has a couple 
of long and short names as attributes:
- long_abbrev
- short_abbrev (one letter abbreviation)
- original_base
- full_name
- modified (True/False)
        """
    parse = PDBParser()
    superimpose = ModernaSuperimposer()

    
    def get_remove_rule(self):
        """
Prepares a rule for removing modification from a residue.
Returns a dict with a rule:
{
'fixed': [link atom names from an original residue],
'moved': [link atom names from a new (standard) base]
}
        """
        if self.long_abbrev in ['Y', 'Ym', 'm3Y', 'm1Y', 'm1acp3Y']: 
            return {'fixed':['C5', 'C4', 'C6'], 'moved':['N1', 'C2', 'C6']}
        elif self.purine: 
            return {'fixed':['N9', 'C8', 'C4'], 'moved':['N9', 'C8', 'C4']}
        elif self.pyrimidine: 
            return {'fixed':['N1', 'C2', 'C6'], 'moved':['N1', 'C2', 'C6']}
        else: 
            raise ModernaResidueError('Residue %d: could not get a removing rule.' % self.number)


    def remove_modification(self):
        """
        Removes a modification from this residue.
        It removes all unnecessary atoms and adds a standard base,
        corresponding to the originating base of the modified one.
        acording to removing modification rule (prepared by get_remove_rule()).
        """
        if not self.modified: 
            raise ModernaResidueError('Residue %d: the residue does not have any modification. Could not remove modification.' %self.number)
        elif self.long_abbrev in ['X', 'Xm']: 
            raise ModernaResidueError('Residue %d: unidentified residue. Could not remove modification.' % self.number)
        elif self.long_abbrev in ['dA', 'dG', 'dC', 'dT']:
            rule = {}
            rule['modification_name'] = self.long_abbrev
            rule['original_base'] = self.original_base
            rule['remove'] = ''
            if self.long_abbrev == 'dT':
                if self.child_dict.has_key('C7'):
                    rule['remove'] = 'C7'
                elif self.child_dict.has_key('C5M'):
                    rule['remove'] = 'C5M'
            rule['moved_link_atoms'] = ["C3'", "C2'", "C1'"]
            rule['fixed_link_atoms'] = ["C3'", "C2'", "C1'"]
            rule['fragment_file_name'] = 'C1pC2pC3p_O2p.pdb'
            rule['pdb_abbrev'] = 'D' + self.long_abbrev[1]
            self.add_single_fragment(rule)
        else:
            new_base = self.parse.get_structure(self.original_base, BASE_PATH + self.original_base + '.ent')[0]['C'][(' ', 54, ' ')]
            triplet_names = self.get_remove_rule()
            
            self.superimpose.get_atoms([self], triplet_names['fixed'], 'fixed')
            self.superimpose.get_atoms([new_base], triplet_names['moved'], 'moved')
            self.superimpose.moved_atoms = new_base.child_list
            self.superimpose.superimpose()

            try:
                for atom in self.child_list[:]:
                    if atom.id not in BACKBONE_RIBOSE_ATOMS: 
                        self.detach_child(atom.id)
                for atom in new_base:
                    self.add(atom)
            except:
                raise ModernaResidueError('Residue %d: couldn not remove unnecesary and add proper atoms' % self.number)

        self.change_name(self.original_base)
        self.modified = False
        self.set_bfactor(B_FACTOR_REMOVE_MODIF)
            

    def get_modification_rules(self, modification_name, separator=' | '):
        """
Prepares a rule for adding a modification. 
Rules describe which fragments add and how to do this
to obtain a residue with given modification.
Returns list of dicts with rules for adding a single fragment.

Keys in each rule dict: ['modification_name', 'original_base', 'remove', 
'moved_link_atoms', 'fixed_link_atoms', 'fragment_file_name', 'pdb_abbrev']
        """
        rules=[]
        # KR: Ooops, is the RULES file really parsed each time when a modification
        #     is added or removed? Time for optimization.
        try: 
            infile = open(ADDING_MODIFICATION_RULES_PATH)
        except IOError:
            raise ModernaResidueError('File does not exist: %s ' % ADDING_MODIFICATION_RULES_PATH)              

        for line in infile:
            line = line.strip().split(separator) 
            if line[0].strip() == modification_name:
                rule = {}
                rule['modification_name'] = line[0]
                rule['original_base'] = line[1]
                rule['remove'] = line[2]
                rule['moved_link_atoms'] = line[3].split(',')
                rule['fixed_link_atoms'] = line[4].split(',')
                rule['fragment_file_name'] = line[5]
                rule['pdb_abbrev'] = line[6]
                rules.append(rule)
        return rules


    def add_single_fragment(self, rule):
        """
        Adds a fragment to a residue.

        Arguments:
        - an adding rule dict (prepared by get_modification_rules())
        """
        try:
            fragment = self.parse.get_structure('fragment', MODIFICATION_FRAGMENTS_PATH + rule['fragment_file_name'])[0]['A'][('H_UNK', 1, ' ')]
        except IOError:
            raise ModernaResidueError('File does not exist: %s' % MODIFICATION_FRAGMENTS_PATH)       

        self.superimpose.get_atoms([self], rule['fixed_link_atoms'], 'fixed')
        self.superimpose.get_atoms([fragment], rule['moved_link_atoms'], 'moved')
        self.superimpose.moved_atoms = fragment.child_list
        self.superimpose.superimpose()
        # KR: This code block also occurs in remove_modification.
        # It could be made a method of superimposer.

        if rule['remove']: 
            delete_atoms = rule['remove'].split(',') + rule['moved_link_atoms']
        else:
            delete_atoms = rule['moved_link_atoms']
        try:
            for atom_name in delete_atoms:
                self.detach_child(atom_name)
            for atom in fragment:
                self.add(atom)       
        except ModernaError, e:
            raise e
        except:
            raise ModernaResidueError('Residue %d: could not remove unnecessary and add proper atoms' % self.number)
            #TODO: remove except:


    def add_modification(self, modification_name):
        """
        Adds a modification to a residue.
        It adds single fragments (add_single_fragment) 
        according to adding modification rules (get_modification_rules).

        Arguments:
        - modification name (as a long abbreviation)
        """
        try:
            if modification_name in [ANY_RESIDUE, MISSING_RESIDUE]: 
                raise ModernaResidueError('Residue %d: expected a modification name, instead got missing/any residue abbreviation "%s"'\
                                % (self.number, modification_name))
            else:
                if self.long_abbrev == UNKNOWN_RESIDUE_SHORT:
                    self.mutate_unknown_residue()
                if self.modified: 
                    self.remove_modification()            
                rules = self.get_modification_rules(modification_name)
                if not rules:
                    raise ModernaResidueError('Residue %d: there is no rule for adding this modification. Check modification name "%s".' \
                        %(self.number, modification_name))
                else:
                    if rules[0]['original_base'] != self.original_base: 
                        self.exchange_base(rules[0]['original_base'])
                    for rule in rules: 
                        self.add_single_fragment(rule)
                    self.change_name(modification_name)   
                    self.set_bfactor(B_FACTOR_ADD_MODIF)                 
        except IOError: 
            raise ModernaResidueError('Residue %d: could not add modification.' % self.number)
       
    def get_exchange_rule(self, new_name):
        """
        Prepares a rule for exchanging a base.
        Returns a dict with the rule.
        {
        'fixed':[link atom names from original residue],
        'moved':[link atom names from new base],
        'remove':[atom names that must be removed from the original residue (old base atoms)]
        }
        """
        rule = {}
        if self.purine: 
            rule['fixed'] = ['N9', 'C4', 'C8']
            rule['remove'] = ['N9', 'C8', 'C4', 'N1', 'C6', 'C5', 'N7', 'C8', 'N3', 'O6', 'N2', 'N6', 'C2']
        elif self.pyrimidine:
            rule['fixed'] = ['N1', 'C2', 'C6']
            rule['remove'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'N4', 'C5', 'C6']
        else: 
            raise ModernaResidueError('Residue %d: could not get exchange rule for name %s' % self.number, new_name)

        if new_name in ['A', 'G']:
            rule['moved'] = ['N9', 'C4', 'C8']
        elif new_name in ['C', 'U']:
            rule['moved'] = ['N1', 'C2', 'C6']
            #TODO: make constants
        else: 
            raise ModernaResidueError('Residue %d: could not get exchange rule for name %s' % self.number, new_name) 
        return rule


    def exchange_base(self, new_name):
        """
        Exchanges standard bases in a residue.

        Arguments:
        - a new base name ('A','G','C' or 'U') 
        """
        if self.long_abbrev == UNKNOWN_RESIDUE_SHORT:
            self.mutate_unknown_residue()
        if self.modified: 
            self.remove_modification() 
        rule = self.get_exchange_rule(new_name)
        new_base = self.parse.get_structure(self.original_base, BASE_PATH+new_name+'.ent')[0]['C'][(' ', 54, ' ')]
            
        self.superimpose.get_atoms([self], rule['fixed'], 'fixed')
        self.superimpose.get_atoms([new_base], rule['moved'], 'moved')
        self.superimpose.moved_atoms = new_base.child_list
        self.superimpose.superimpose()
        
        for atom in rule['remove']:
            if atom in self.child_dict.keys(): 
                self.detach_child(atom)
        for atom in new_base: 
            self.add(atom) 
        self.change_name(new_name)
        self.set_bfactor(B_FACTOR_EXCHANGE)

    def get_renumbered_resi(self, number):
        resi = ModernaResidue(self, alphabet_entry=self.alphabet_entry)
        resi.change_number(number)
        return resi

    def make_backbone_only_residue(self, include_N=True):
        """
        Cuts the base out of a residue so only backbone and ribose atoms stay.
        """
        backbone_atoms = BACKBONE_RIBOSE_ATOMS[:]
        if include_N: 
            nglyc = self["N*"]
            backbone_atoms += [nglyc.name.strip()]
        for atom in self.child_list[:]:
            # copying list because list is modified while iterated
            if atom.name.strip() not in backbone_atoms: 
                self.detach_child(atom.id)

        for at_name in backbone_atoms:
            if at_name not in self.child_dict.keys():
                raise ModernaResidueError('Residue %s: backbone is not complete. Missing atom %s' %(self.number, at_name))
        self.change_name(ANY_RESIDUE)


    def mutate_unknown_residue(self):
        """
        Makes a C out of unknown residue (X, .) on ribose and N* (N9,N1) atoms.
        When ribose and N* are nor present raises an error.
        C can be then changed into any modification / standard residues.
        """
        for x in RIBOSE:
            if not self.has_id(x): 
                raise ModernaResidueError('Residue %s: cannot mutate unknown residue. Missing ribose atom %s' %(self.number, x))
        if not self.has_id('N9') and not self.has_id('N1'):
            raise ModernaResidueError('Residue %s: cannot mutate unknown residue. Missing ribose atom N* (N1 or N9)' %(self.number))
        
        try:
            fragment = self.parse.get_structure('fragment', MODIFICATION_FRAGMENTS_PATH +'ribose_C.pdb')[0]['A'][(' ', 1, ' ')]
        except IOError:
            raise ModernaResidueError('File does not exist: %s' % MODIFICATION_FRAGMENTS_PATH+'ribose_C.pdb')       

        self.superimpose.get_atoms([self], ["O4'", "C1'", "C2'"], 'fixed')
        self.superimpose.fixed.append(self["N*"])
        self.superimpose.get_atoms([fragment], ["O4'", "C1'", "C2'", 'N1'], 'moved')
        self.superimpose.moved_atoms = fragment.child_list
        self.superimpose.superimpose()
        
        for x in list(self): # must copy items before deleting.
            if  x.name not in BACKBONE_RIBOSE_ATOMS:
               self.detach_child(x.name)
               
        for x in fragment:
            if x.name not in  ["O4'", "C1'", "C2'"]:
                self.add(x)
                
        self.change_name('C')   
        
        
    def mutate(self, abbr):
        """
        Deals with exchanging base, adding and removing modification.
        Decides by itself which operation to perform.

        Arguments:
        - a long abreviation ('A', 'mA')
        """
        if abbr == UNKNOWN_RESIDUE_SHORT: 
            self.make_backbone_only_residue()
        if self.long_abbrev != abbr:
            if self.long_abbrev == UNKNOWN_RESIDUE_SHORT:
               self.mutate_unknown_residue()
            if self.modified: 
                self.remove_modification()
            if abbr not in ['A', 'U', 'G', 'C']: 
                self.add_modification(abbr)
            else: 
                self.exchange_base(abbr)

    def set_bfactor(self, b_value):
        """Sets the same b factor for all atoms in the residue"""
        for atom in self:
            atom.set_bfactor(b_value)

    def rotate_chi(self, angle=90):
        """
        Rotates the base around the glycosidic bond
        by the given chi angle
        """
        #TODO: is this method used?
        try:
            vC1 = self["C1'"].get_vector()
            vN = self["N*"].get_vector()
        except KeyError:
            raise ModernaResidueError("Residue %d: could not find atom C1' or N9 (if A or G)or N1 (if C or U)." % self.number ) 

        v_glycosidic = vN - vC1
        matrix = rotaxis2m(math.radians(angle), v_glycosidic)

        for atom in self:
            if atom.name not in BACKBONE_RIBOSE_ATOMS and atom != self["N*"]:
                v_atom = self[atom.name].get_vector()
                vec = v_atom - vC1
                vec = vec.left_multiply(matrix)
                vec += vC1
                atom.set_coord([vec[0], vec[1], vec[2]])
                
    #
    # helper for SimRNA
    #
    def lock(self):
        """Sets all occupancies to zero."""
        for atom in self.child_list:
            atom.occupancy = 0.0
    
    def unlock(self):
        """Sets all occupancies to one."""
        for atom in self.child_list:
            atom.occupancy = 1.0
            
    def b_zero(self):
        """Sets B-factor to zero"""
        for atom in self.child_list:
            atom.bfactor = 0.0
        


#
# AUXILIARY METHOD FOR TESTING PURPOSES
#
#TODO: move somewhere else.
def write_hydrogens(hydrogens, filename='hydro.pdb'):
    """Writes the hydrogen as a PDB file"""
    from Bio.PDB import PDBIO
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure
    import numpy as np
    my_structure = Structure('my structure id')
    my_model = Model(0)
    my_chain = Chain('A')
    my_structure.add(my_model)
    my_model.add(my_chain)
    # add atoms
    for i, coord in enumerate(hydrogens):
        my_residue  = Residue((' ', i, ' '), 'UNK', ' ')
        coord_array = np.array(list(coord),'f')
        atom = Atom('H', coord_array, 0.0, 1.0, ' ', 'H', i+1, element='H')
        my_residue.add(atom)
        my_chain.add(my_residue)
    # Save structure.
    out = PDBIO()
    out.set_structure(my_structure)
    out.save(filename)
        
