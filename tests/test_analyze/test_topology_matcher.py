#!/usr/bin/env python
#
# test_topology_matcher.py
#
# Tests the topology matching algorithm.
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

from unittest import main, TestCase
from moderna.analyze.MolecularGraph import AnnotatedMolecule
from moderna.analyze.MolGraphParser import Molecule
from test_data import TEST_DATA_PATH
import os

ATTRIBUTES = ['h-bond acceptor','h-bond donor','h-bond hydrogen',\
            'sp2-hybrid','sp3-hybrid','chiral_atom','conjugated_system'
            ]

class MolParserTests(TestCase):

    def test_init(self):
        """A new molecule should have no atoms."""
        m = Molecule()
        self.assertEqual(len(m),0)
        
    def check_molecule(self, m, elements):
        """Checks element numbers and makes sure each atom has basic attributes."""
        count = {}
        for atom in m:
            count.setdefault(atom.element,0)
            count[atom.element] += 1
        for key in elements:
            self.assertEqual(count.get(key,0),elements[key])
        for atom in m:
            self.assertTrue(atom.has_key('coordinates'))
            self.assertTrue(atom.has_key('atom_name'))
            self.assertTrue(atom.has_key('attributes'))
            self.assertTrue(atom.has_key('rdf_index'))    
        
    def test_parse_molfile(self):
        """Resulting molecules should have the right number of atoms for each element."""
        # water molecule from kegg
        m = Molecule()
        m.parse_molfile(TEST_DATA_PATH+'topology_data/c00001.mol')
        self.assertEqual(len(m),3)
        self.check_molecule(m, {'H':2,'O':1})
        # ATP molecule from kegg
        m = Molecule()
        m.parse_molfile(TEST_DATA_PATH+'topology_data/c00002.mol')
        self.assertEqual(len(m),47)
        self.check_molecule(m, {'N':5,'C':10,'O':13,'P':3,'H':16})
        # Biopath molecule
        m = Molecule()
        m.parse_molfile(TEST_DATA_PATH+'topology_data/biopath_single.mol')
        self.assertEqual(len(m),24)
        self.check_molecule(m, {'C':5,'O':3,'N':3,'H':13})
        # Biopath molecule
        m = Molecule()
        m.parse_molfile(TEST_DATA_PATH+'topology_data/biopath_double.mol')
        self.assertEqual(len(m),12)
        self.check_molecule(m, {'C':3,'O':3,'H':6})
    
    def test_bonds(self):
        """bonds should be in place"""
        # water molecule from kegg
        m = Molecule()
        m.parse_molfile(TEST_DATA_PATH+'topology_data/c00001.mol')
        h = []
        o = None
        for atom in m:
            if atom.element == 'H': h.append(atom)
            else: o = atom
        self.assertEqual(len(o.bonds),2)
        self.assertEqual(len(h[0].bonds),1)
        self.assertEqual(len(h[1].bonds),1)
        for bond in o.bonds:
            self.assertTrue(bond.atom2 in h)
            self.assertEqual(bond.valence,1)
        self.assertEqual(h[0].bonds[0].atom2,o)
        self.assertEqual(h[0].bonds[0].valence,1)
        self.assertEqual(h[1].bonds[0].atom2,o)
        self.assertEqual(h[1].bonds[0].valence,1)
    
    def test_bondschemes(self):
        """Check bonds using bonding scheme calculation."""
        m = Molecule()
        m.parse_molfile(TEST_DATA_PATH+'topology_data/biopath_double.mol')
        for atom in m:
            atom.detect_bonding_schema()
        bs = [atom['bondschema'] for atom in m]
        bondschemes = ['CC1H111','CC11H1O1','CC1O12','OC1H1','OC2']
        for b in bondschemes:
            self.assertTrue(b in bs)
        
        
class AnnotatedMoleculeTests(TestCase):

    def check_molecule(self,m):
        """
        Makes sure the attributes for each atom are set.
        """
        for atom in m:
            self.assertTrue(atom.has_key('oxidative_number'))
            self.assertTrue(atom.has_key('partial_charge'))
            self.assertTrue(atom.has_key('protortype'))
            self.assertTrue(atom.has_key('wang_atomtype'))
            self.assertTrue(atom.has_key('charge'))
            self.assertTrue(atom.has_key('free_electron_pairs'))
            self.assertTrue(atom.has_key('electronegativity'))
            
        
    def test_attributes(self):
        """Chemical attributes should be all set."""
        m = AnnotatedMolecule()
        m.parse_molfile(TEST_DATA_PATH+'topology_data/biopath_double.mol')
        m.annotate()
        for a in ATTRIBUTES:
            exists = False
            for atom in m:
                # print atom['attributes']
                if a in atom['attributes']:                    
                    exists = True
            self.assertTrue(exists)
            
    def test_chemical_groups(self):
        """Chemical groups should be recognized correctly."""
        # simple example: water
        groups = {
            'OH11':['o_water'],
            'HO1':['h_water'],
            }
        self.check_groups(TEST_DATA_PATH+'topology_data/c00001.mol', groups)
        # complicated example: lactate
        groups = {
            'CC1H111':['alpha-hydroxyl-c_methyl', 'alpha-secondary_alcohol-c_methyl', 'c_methyl'],
            'CC11H1O1':['alpha-carboxy-c_secondary_alcohol', 'alpha-methyl-c_secondary_alcohol', 'c_secondary_alcohol'],
            'OC1H1':('alt',['o_hydroxyl'],['o_carboxyl']),
            'CC1O12':['alpha-hydroxyl-c_carboxyl', 'alpha-secondary_alcohol-c_carboxyl', 'c_carboxyl'],
            'OC2':['o_carboxyl'],
            'HC1':['h_other_hydrogen'],
            'HO1':['h_other_hydrogen'],
            }
        self.check_groups(TEST_DATA_PATH+'topology_data/biopath_double.mol', groups)
        
    def check_groups(self, filename, groups):
        m = AnnotatedMolecule()
        m.parse_molfile(filename)
        m.annotate()
        for atom in m:
            at = atom['attributes']
            group = groups[atom['bondschema']]
            for an in ATTRIBUTES: 
                if an in at: at.remove(an)
            at.sort()
            if group[0] == 'alt':
                for g in group[1:]:
                    if g == at: 
                        group = g
                        break
            self.assertEqual(at,group)
        
    '''
    This test might be useful in Modomics.
    It checks whether all .mol files are recognized correctly.
    (taken out of ModeRNA, because BaseRecognizer tests this now.).
    KR (2011/04/01)
    
    def test_detect_modified_bases(self):
        """All structures from Modomics should be recognized correctly."""
        ok = 0
        total = 0
        for fn in os.listdir(TEST_DATA_PATH+'topology_data/modified_bases/'):
            if fn[-4:] == '.mol':
                m = AnnotatedMolecule()
                m.parse_molfile(TEST_DATA_PATH+'topology_data/modified_bases/'+fn)
                if len(m)==0: 
                    # skip the dummy molfiles
                    continue
                modname = fn[:-4]
                modifications = m.detect_modified_bases()
                if 'pyrimidine' in modifications: modifications.remove('pyrimidine')
                if 'purine' in modifications: modifications.remove('purine')
                if 'riboside' in modifications: modifications.remove('riboside')
                if 'phosphate' in modifications: modifications.remove('phosphate')
                if modname in ['galQtRNA','manQtRNA','T','D','m5D']:
                    ok += 1 # indistinguishable by the algorithm.
                
                elif len(modifications)==1 and modifications[0] == modname:
                    ok += 1
                else:
                    print "modification %s not recognized. check mol file. %s"%(modname, str(modifications))
                total += 1
        print "\nModified bases examined   : %i"%total
        print "Modified bases recognized : %i"%ok
        
        self.assertEqual(ok, total)
        '''
    

if __name__== '__main__':
    main()
