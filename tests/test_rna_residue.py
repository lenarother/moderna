#!/usr/bin/env python
#
# test_copy_residue.py
#
# unit tests for residue duplication functionality
#
__author__ = "Magdalena Musielak, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Musielak"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Prototype"

from unittest import main, TestCase
from moderna.RNAResidue import RNAResidue
from moderna.analyze.BaseRecognizer import BaseRecognizer, BaseRecognitionError
from Bio.PDB import PDBParser
from moderna.util.Errors import RNAResidueError
from moderna.sequence.ModernaAlphabet import Alphabet
from moderna.Constants import BIO153

from test_data import *

class RNAResidueTests(TestCase):
    def setUp(self):
        """Loads the A residue to start with."""
        self.a=PDBParser().get_structure('test_struc',A_RESIDUE)[0].child_list[0].child_list[0]
        self.chain=PDBParser().get_structure('test_struc',MINI_TEMPLATE)[0].child_list[0]
        self.chain2=PDBParser().get_structure('test_struc',MINI_TEMPLATE)[0].child_list[0]
    
    def tearDown(self):
        self.a = None
        self.chain = None
        self.chain2 = None

    def test_residue_identity(self):
        """Moderna residues need to be discinct by __eq__ unless they are the same object."""
        r1 = RNAResidue(self.chain.child_list[2])
        r2 = RNAResidue(self.chain2.child_list[2])
        r3 = RNAResidue(self.chain.child_list[3])
        r4 = RNAResidue(self.chain.child_list[2])
        self.assertEqual(r1, r1)
        self.assertNotEqual(r1, r2)
        self.assertNotEqual(r1, r3)
        self.assertNotEqual(r1, r4)
        
    def test_atom_parent(self):
        """Atoms should link back to their parent."""
        resi = RNAResidue(self.a)
        for atom in resi:
            self.assertEqual(atom.get_parent(),resi)
            
    def test_atom_name(self):
        """Atoms should have the right names."""
        resi = RNAResidue(self.a)
        self.assertEqual(resi["C4'"].name,"C4'")
        self.assertEqual(resi["C4'"].fullname," C4'")
        if BIO153:
            self.assertEqual(resi["C4'"].element,'C')
            
    def test_init(self):
        """Residues should be initializable."""
        a = RNAResidue(self.a)
        self.assertEqual(a.identifier,'1')
        a = RNAResidue(self.chain.child_list[0])
        self.assertEqual(a.identifier,'1')
        a = RNAResidue(self.chain.child_list[-1])
        self.assertEqual(a.identifier,'15')
        
    def test_init_recog_base(self):
        """Base recognition should succeed regardless of parameter"""
        alphabet = Alphabet()
        # recognition with base recognizer
        a = RNAResidue(self.chain.child_list[9])
        self.assertEqual(a.long_abbrev,'m2G')
        # assignment by alphabet entry
        # IMPORTANT FOR BASE RECOGNIZER BYPASS
        a = RNAResidue(self.chain.child_list[9], alphabet['m2G'])
        self.assertEqual(a.long_abbrev,'m2G')
        # assignment by wrong alphabet entry should succeed!
        a = RNAResidue(self.chain.child_list[9], alphabet['m7G'])
        self.assertEqual(a.long_abbrev,'m7G')
        
    def test_renumber(self):
        res = RNAResidue(self.chain.child_list[2])
        res.change_number('64')
        self.assertEqual(res.identifier, '64')
        self.assertEqual(res.id, (' ',64, ' '))
        res.change_number('133E')
        self.assertEqual(res.identifier, '133E')
        self.assertEqual(res.id, (' ',133, 'E'))
        res.change_number('2')
        self.assertEqual(res.identifier, '2')
        self.assertEqual(res.id, (' ',2, ' '))
        
    def test_glycosidic_n(self):
        """Finds N* in tough cases."""
        chain = PDBParser().get_structure('test_struc', 'test_data/gaps/1h3e_B.pdb')[0].child_list[0]
        resi1 = RNAResidue(chain[(' ', 15, ' ')])
        self.assertTrue(resi1['N*'])
        resi2 = RNAResidue(chain[(' ', 16, ' ')])
        self.assertRaises(RNAResidueError, resi2.__getitem__, 'N*')

    def test_purine(self):
        """RNAResidue recognizes purines."""
        res = RNAResidue(self.chain.child_list[0])
        self.assertTrue(res.purine)
        res =RNAResidue(self.chain.child_list[9])
        self.assertTrue(res.purine)
        res =RNAResidue(self.chain.child_list[1])
        self.assertFalse(res.purine)
        res = RNAResidue(self.chain.child_list[11])
        self.assertFalse(res.purine)
        
    def test_pyrimidine(self):
        """ModernaResidue recognizes pyrimidines."""
        res = RNAResidue(self.chain.child_list[1])
        self.assertTrue(res.pyrimidine)
        res = RNAResidue(self.chain.child_list[11])
        self.assertTrue(res.pyrimidine)
        res = RNAResidue(self.chain.child_list[0])
        self.assertFalse(res.pyrimidine)
        res = RNAResidue(self.chain.child_list[9])
        self.assertFalse(res.pyrimidine)

if __name__ == '__main__':
    main()
  
