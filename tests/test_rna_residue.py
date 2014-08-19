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
        
    # -------- TESTS FOR INTEGRITY ----------------
    def test_is_backbone_complete(self):
        """Missing backbone atoms are recognized."""
        for resi in self.chain:
            resi = RNAResidue(resi)
            #self.assertTrue(resi.is_backbone_complete())
        # .. and now the negative example.
        chain=PDBParser().get_structure('test_struc',INCOMPLETE_BACKBONE)[0].child_list[0]
        for resi in chain:
            resi = RNAResidue(resi)
            self.assertFalse(resi.is_backbone_complete())
        
    def test_is_backbone_intact(self):
        """Should check all kinds of backbone discontinuities in one residue."""
        chain=PDBParser().get_structure('test_struc',BROKEN_BACKBONE)[0].child_list[0]
        residues = [r for r in chain]
        for resi in residues[:5]:
            mr = RNAResidue(resi)
            self.assertFalse(mr.is_backbone_intact())
        mr = RNAResidue(chain[6])
        self.assertTrue(mr.is_backbone_intact())
        
    def test_is_backbone_intact_5p3p(self):
        """Should check all kinds of backbone discontinuities in one residue."""
        chain=PDBParser().get_structure('test_struc',BROKEN_BACKBONE)[0].child_list[0]
        residues = [r for r in chain]
        result_5p = []
        result_3p = []
        for resi in residues[:6]:
            mr = RNAResidue(resi)
            result_5p.append(mr.is_backbone_intact(mode="5'"))
            result_3p.append(mr.is_backbone_intact(mode="3'"))
        self.assertEqual(result_5p, [False, False, False, False, True, True])
        self.assertEqual(result_3p, [True, True, True, False, False, True])
        
    def test_is_phosphate_intact(self):
        """Should check whether OP1 and OP2 are in place"""
        chain=PDBParser().get_structure('test_struc','test_data/rna_structures/bb_messed_up.pdb')[0].child_list[0]
        resi1 = RNAResidue(chain[('H_c  ', 32, ' ')])
        resi2 = RNAResidue(chain[(' ', 33, ' ')])
        self.assertTrue(resi1.is_phosphate_intact())
        self.assertFalse(resi2.is_phosphate_intact())
    
    def test_is_backbone_congested(self):
        """Should check whether backbone atoms clash into rest of the structure."""
        resi = RNAResidue(self.chain.child_list[2])
        self.assertFalse(resi.is_backbone_congested())
        # now check a structure where the backbone clashes into O2'
        chain=PDBParser().get_structure('test_struc','test_data/rna_structures/bb_messed_up.pdb')[0].child_list[0]
        resi = RNAResidue(chain[('H_c  ', 32, ' ')])
        self.assertTrue(resi.is_backbone_congested())
        resi = RNAResidue(chain[(' ', 33, ' ')])
        self.assertTrue(resi.is_backbone_congested())


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
  
