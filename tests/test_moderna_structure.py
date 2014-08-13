#!/usr/bin/env python
#
# test_moderna_structure.py
#
# unit tests for ModernaStructure class
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

from unittest import main, TestCase
from moderna.ModernaStructure import ModernaStructure, ModernaResidue
import os, tempfile
from Bio.PDB.Model import Model
from moderna.Errors import ModernaStructureError, ModernaResidueError, RNAChainError
from moderna.ModernaSequence import Sequence
from moderna.LogFile import log
from moderna import load_model

from test_data import *

OUTPUT = 'test_data/test_output.ent'

class ModernaStructureTests(TestCase):
    
    def tearDown(self):
        if os.path.exists('temp.pdb'):
            os.remove('temp.pdb')

    def test_add_write(self):
        """A modified base should be written to PDB."""
        s=ModernaStructure('file',A_RESIDUE)
        s['1'].add_modification('m1A')
        s.write_pdb_file(OUTPUT)
        t = ModernaStructure('file', OUTPUT)
        self.assertEqual(t['1'].long_abbrev,s['1'].long_abbrev)
        
    def test_find_residues_in_range(self):
        """Should return a list of residues."""
        s = ModernaStructure('file',MINI_TEMPLATE)
        result = s.find_residues_in_range('3', '9')
        self.assertEqual(result, ['4', '5', '6', '7', '8'])
        
    def test_renumber_residue(self):
        s = ModernaStructure('file', MINI_TEMPLATE)
        self.assertTrue(s['3'])
        self.assertRaises(RNAChainError,s.__getitem__,'30')
        s.renumber_residue('3','30')
        self.assertTrue(s['30'])
        self.assertRaises(RNAChainError,s.__getitem__,'3')
        self.assertEqual(s['30'].identifier,'30')

    def test_get_renumbered_resi(self):
        """Should renumber but not insert a residue."""
        s = ModernaStructure('file', MINI_TEMPLATE)
        resi = s['3'].get_renumbered_resi('100B')
        self.assertNotEqual(resi, s['3'])
        self.assertEqual(resi.identifier, '100B')
        self.assertEqual(s['3'].identifier, '3')

    def test_renumber_residue_letter(self):
        s = ModernaStructure('file', MINI_TEMPLATE)
        self.assertTrue(s['3'])
        self.assertRaises(RNAChainError,s.__getitem__,'3A')
        s.renumber_residue('3','3A')
        self.assertTrue(s['3A'])
        self.assertRaises(RNAChainError,s.__getitem__,'3')
        self.assertEqual(s['3A'].identifier,'3A')

    def test_renumber_residue_with_error(self):
        """Renumbering errors should not remove a residue."""
        s = ModernaStructure('file',MINI_TEMPLATE)
        self.assertRaises(RNAChainError,s.renumber_residue,'5','a')
        s.renumber_residue('5','12A')
        self.assertTrue(s['12A'])
                    
    def test_chain_renumbering(self):
        """Chain should be renumbered"""
        m = ModernaStructure('file',RNA_1EHZ)
        m.renumber_chain('100')
        # check residues
        self.assertEqual(m['100'].long_abbrev,'G')
        self.assertEqual(m['107'].long_abbrev,'U')
        self.assertEqual(m['109'].long_abbrev,'m2G')
        self.assertEqual(m['175'].long_abbrev,'A')
        self.assertEqual(m['115'].long_abbrev,'D')
        self.assertRaises(RNAChainError,m.renumber_chain,'77A')  

    def test_checking_letters_in_residues_numeration(self):
        """Presence of letters in residues numeration should be detected"""
        m = ModernaStructure('file',RNA_1EHZ)
        self.assertFalse(m.check_letters_in_residue_numeration())

    def test_fix_backbone(self):
        s = ModernaStructure('file', FIXABLE_BACKBONE)
        self.assertEqual(s.get_sequence(), Sequence('ACUG_UG'))
        s.fix_backbone()
        self.assertEqual(s.get_sequence(), Sequence('ACUGUG'))
        
    def test_get_modified_resis(self):
        s = ModernaStructure('file', RNA_1EHZ, 'A')
        resis = s.get_modified_residues()
        self.assertEqual(len(resis), 14)
        # water and ions should not be recognized
        s = ModernaStructure('file', RNA_1EHZ, ' ')
        resis = s.get_modified_residues()
        for r in resis:
            print r, s[r]
        self.assertEqual(len(resis), 0)
        
    def test_read_write_phosphate(self):
        """extra 2'3' phosphate is recognized and read correctly"""
        struc=ModernaStructure('file', THREEPRIME_PHOSPHATE, 'B')
        struc.write_pdb_file('temp.pdb')
        struc = ModernaStructure('file', 'temp.pdb', 'B')
        self.assertEqual(struc['190'].long_abbrev, '23pA')

    def test_change_sequence(self):
        """Checks whether sequence is changed into given by user"""
        m = ModernaStructure('file',MINI_TEMPLATE) #residues_mixed.ent')
        m.change_sequence('UUUUUUUUUUUUUUU')
        self.assertEqual(m.get_sequence(), Sequence('UUUUUUUUUUUUUUU'))
        self.assertRaises(ModernaStructureError,m.change_sequence,'AAAAA')
        m.change_sequence('GGGG', '1','4')
        self.assertEqual(m.get_sequence(), Sequence('GGGGUUUUUUUUUUU'))
        m.change_sequence('AAAA', '12')
        self.assertEqual(m.get_sequence(), Sequence('GGGGUUUUUUUAAAA'))
        m.change_sequence('CC', None,'2')
        self.assertEqual(m.get_sequence(), Sequence('CCGGUUUUUUUAAAA'))
        m.change_sequence(Sequence('UUUUUUUUUUUUUUU'))
        self.assertEqual(m.get_sequence(), Sequence('UUUUUUUUUUUUUUU'))
        unk = load_model(PDB_UNK)
        unk.change_sequence('DAAAPAEA?A7A')
        self.assertEqual(unk.get_sequence(), Sequence('DAAAPAEA?A7A'))
        
    def _test_messy_atoms(self):
        """Structure with messed up atoms is rejected."""
        #TEST DISABLED
        before = len(log.contents)
        m = ModernaStructure('file',MESSED_ATOMS) 
        after = len(log.contents)
        self.assertTrue(after>before)
        self.assertTrue(re.find('no nucleotide residues found', log.contents[-1]))
        #TODO: new test, code not implemented yet.
        
    def test_cap5(self):
        """Caps in mRNA are recognized."""
        struc = ModernaStructure('file', CAP_EXAMPLE, 'B')
        seq = str(struc.get_sequence())
        self.assertEqual(seq, 'x_GAA')
        self.assertEqual(struc['900'].long_abbrev, 'm7Gpp_cap')

if __name__ == '__main__':
    main()
    
