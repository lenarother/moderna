#!/usr/bin/env python
"""
Unit Tests for RNAChain class
"""

from unittest import main, TestCase
from moderna.RNAResidue import RNAResidue
from moderna.RNAChain import RNAChain
from Bio.PDB.PDBParser import PDBParser
from moderna.util.Errors import RNAChainError
from moderna.sequence.ModernaSequence import Sequence

from moderna.tests.test_data import *

OUTPUT = 'test_output.ent'


class RNAChainTests(TestCase):

    def test_init_empty(self):
        """Empty chain can be initialized."""
        s = RNAChain()

    def test_init_struc(self):
        struc = PDBParser().get_structure('test', MINI_TEMPLATE)
        s = RNAChain('structure', struc, seq=Sequence("GCGGAUUUALCUCAG"))
        self.assertEqual(s.get_sequence(), Sequence("GCGGAUUUALCUCAG"))

    def test_init_chain(self):
        struc = PDBParser().get_structure('test', MINI_TEMPLATE)
        s = RNAChain('chain', struc[0]['A'], seq=Sequence("GCGGAUUUALCUCAG"))
        self.assertEqual(s.get_sequence(), Sequence("GCGGAUUUALCUCAG"))

    def test_init_residues(self):
        struc = PDBParser().get_structure('test', MINI_TEMPLATE)
        s = RNAChain('residues', struc[0]['A'].child_list, seq=Sequence("GCGGAUUUALCUCAG"))
        self.assertEqual(s.get_sequence(), Sequence("GCGGAUUUALCUCAG"))

    def test_init_with_sequence(self):
        """Allows to skip sequence recognition"""
        s = RNAChain('file', MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
        self.assertEqual(s.get_sequence(),Sequence("GCGGAUUUALCUCAG"))
        # wrong sequence should also be applied
        #TODO: check whether this is OK
        s = RNAChain('file', MINI_TEMPLATE, seq=Sequence("UUULLYYDDAAGAGA"))
        self.assertEqual(s.get_sequence(),Sequence("UUULLYYDDAAGAGA"))

    def test_true(self):
        """Should always be true."""
        s = RNAChain()
        self.assertTrue(s)

    def test_load_iter(self):
        """Should load PDB files and iterate over all residues."""
        s = RNAChain('file',MINI_TEMPLATE)
        i = 0
        for resi in s:
            i += 1
        self.assertEqual(i,15)

    def test_load_error(self):
        """Should raise exception if the file does not exist."""
        self.assertRaises(RNAChainError, RNAChain, 'file','nonexisting_file.ent')

    def test_load_chain_error(self):
        """Should raise exception if the chain does not exist."""
        self.assertRaises(RNAChainError, RNAChain, 'file', MINI_TEMPLATE, 'B')

    def test_iterate(self):
        """Iteration over RNAResidues with correct numbers."""
        s = RNAChain('file', MINI_TEMPLATE)
        i = 1
        for resi in s:
            self.assertTrue(isinstance(resi, RNAResidue))
            self.assertEqual(resi.id[1], i)
            i += 1

    def test_write_pdb(self):
        """write_pdb should create a file."""
        s = RNAChain('file', MINI_TEMPLATE)
        if os.access(OUTPUT, os.F_OK):
            os.remove(OUTPUT)
        s.write_pdb_file(OUTPUT)
        self.assertTrue(os.access, os.F_OK)

    def test_write_load_sanity(self):
        """After re-reading written file, there should be the same residues."""
        s = RNAChain('file', MINI_TEMPLATE)
        if os.access(OUTPUT, os.F_OK):
            os.remove(OUTPUT)
        s.write_pdb_file(OUTPUT)
        t = RNAChain('file', OUTPUT)
        i = 0
        while i<15:
            i += 1
            self.assertEqual(s[str(i)].id[1],t[str(i)].id[1])
            self.assertEqual(s[str(i)].long_abbrev,t[str(i)].long_abbrev)   

        s_atoms= sum([len(r) for r in s])
        t_atoms= sum([len(r) for r in t])
        self.assertEqual(s_atoms,t_atoms)

    def test_getitem(self):
        """Residues should be accessible by their number as string."""
        s  = RNAChain('file', MINI_TEMPLATE)
        self.assertEqual(s['1'].long_abbrev, 'G')
        self.assertEqual(s['8'].long_abbrev,'U')
        self.assertEqual(s['9'].long_abbrev,'A')
        self.assertEqual(s['12'].long_abbrev,'U')
        self.assertEqual(s['15'].long_abbrev,'G')
        self.assertRaises(RNAChainError,s.__getitem__,'17')

    def test_getitem_number(self):
        """Calling getitem with a number gives an error."""
        s = RNAChain('file', MINI_TEMPLATE)
        self.assertRaises(RNAChainError,s.__getitem__,1)

    def test_getitem_slice(self):
        """Residues should be accessible by slice."""
        s = RNAChain('file', MINI_TEMPLATE)
        struc_slice = s['2':'5']
        self.assertEqual(len(struc_slice),4)
        self.assertEqual(struc_slice[0].identifier,'2')
        self.assertEqual(struc_slice[3].identifier,'5')
    
    def test_getitem_slice_yourself(self):
        """Residues should be accessible by slice."""
        s = RNAChain('file', MINI_TEMPLATE)
        struc_slice = s['2':'2']
        self.assertEqual(len(struc_slice),1)
        self.assertEqual(struc_slice[0].identifier,'2')
        struc_slice = s['14':'14']
        self.assertEqual(len(struc_slice),1)
        self.assertEqual(struc_slice[0].identifier,'14')
        
    def test_getitem_slice_with_letters(self):
        """Residues should be accessible by slice."""
        s = RNAChain('file',BAD_LOOP_INSERTION)
        struc_slice = s['944':'961']
        self.assertEqual(len(struc_slice), 9)
        self.assertRaises(RNAChainError, s.__getitem__, slice('0', '961'))
        
    
    def test_len(self):
        """should return correct number of residues."""
        s = RNAChain('file', MINI_TEMPLATE)
        self.assertEqual(len(s),15)

    def test_get_sequence(self):
        s = RNAChain('file', MINI_TEMPLATE)
        seq = s.get_sequence()
        self.assertEqual(seq,Sequence("GCGGAUUUALCUCAG"))
        
    def test_get_sequence_with_hydromods(self):
        """Modifications with hydrogens should be recognized."""
        s = RNAChain('file', MOD_WITH_HYDROGENS, 'B')
        seq = s.get_sequence()
        self.assertEqual(seq,Sequence("UPA"))
        
    def test_1ehz(self):
        """Everything should work for a more complicated structure."""
        m = RNAChain('file',RNA_1EHZ)
        # atom count
        self.assertEqual(len(m),76)
        # write
        if os.access(OUTPUT,os.F_OK): os.remove(OUTPUT)
        m.write_pdb_file(OUTPUT)
        self.assertTrue(os.access(OUTPUT,os.F_OK))
        # reload
        n = RNAChain('file',OUTPUT)
        self.assertEqual(len(n),76)
        # check iter
        res = [r for r in n]
        self.assertEqual(len(res),76)
        # check residues
        self.assertEqual(n['1'].long_abbrev,'G')
        self.assertEqual(n['8'].long_abbrev,'U')
        self.assertEqual(n['10'].long_abbrev,'m2G')
        self.assertEqual(n['76'].long_abbrev,'A')
        self.assertEqual(n['16'].long_abbrev,'D')
        self.assertRaises(RNAChainError,n.__getitem__,'77')        

    def test_sort_residues(self):
        m = RNAChain('file',MINI_TEMPLATE) #residues_mixed.ent')
        m.sort_residues()
        residues=[]
        for resi in m:
            residues.append(resi.id[1])
        residues_copy=residues[:]
        residues_copy.sort()
        self.assertEqual(residues,residues_copy)
    
    def test_PDB_file_with_negative_numbering(self):
        """Checks whether the order of residues with negative NRs is correct"""
        # ...because this was the major problem with negative numbering...
        # ...that ModeRNA read residues in a wrong order!
        s = RNAChain('file', NEG_NUM_PDB, 'E')
        self.assertEqual(s.get_sequence().seq_with_modifications, \
        'GGCUGCCUGGGUCCGCCUUGAGUGCCCGGGUGAGAAGCAUGAUCCCGGGUAAUUAUGGCGGACCCACA')
        
        self.assertEqual(s.get_sequence().seq_without_modifications, \
        'GGCUGCCUGGGUCCGCCUUGAGUGCCCGGGUGAGAAGCAUGAUCCCGGGUAAUUAUGGCGGACCCACA')

    def test_ignore_empty_residues(self):
        """Check whether residues with only one H are read and written."""
        s = RNAChain('file', H_ONLY, ' ')
        self.assertEqual(len(s), 0)
        
    def test_first_resi(self):
        """Should return firs residue"""
        s = RNAChain('file', MINI_TEMPLATE)
        self.assertEqual(s.first_resi.identifier, '1')
        
    def test_last_resi(self):
        """Should return firs residue"""
        s = RNAChain('file', MINI_TEMPLATE)
        self.assertEqual(s.last_resi.identifier, '15')
        
    def test_get_region(self):
        """Checks whether a region can be extracted from structure"""
        s = RNAChain('file', MINI_TEMPLATE)
        s2 = s.get_region('3', '6')
        self.assertTrue(isinstance(s2, RNAChain))
        self.assertEqual(len(s2.moderna_residues), 4)
        self.assertEqual(s2.get_sequence(), Sequence('GGAU'))
        s3 = s.get_region('3')
        self.assertEqual(s3.get_sequence(), Sequence('GGAUUUALCUCAG'))
        s4 = s.get_region(stop_id='6')
        self.assertEqual(s4.get_sequence(), Sequence('GCGGAU'))


class SpeedTests(TestCase):
    """
    This is not part of the canonical UnitTests. Use it if you want to check performance.
    """
    def test_slow(self):
        for i in range(10):
            m = RNAChain('file',RNA_1EHZ)
    
    def test_fast(self):
        for i in range(10):
            s = Sequence('GCGGAUUUALCUCAGDDGGGAGAGCRCCAGABU#AAYAP?UGGAG7UC?UGUGTPCG"UCCACAGAAUUCGCACCA')
            m = RNAChain('file',RNA_1EHZ, seq=s)
        

if __name__ == '__main__':
    main()
    
