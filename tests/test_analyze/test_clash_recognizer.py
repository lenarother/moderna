#!/usr/bin/env python
#
# test_find_clashes.py
#
# unit tests for find_clashes module
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Tomasz Puton"
__copyright__ = "Copyright 2008, The Moderna Project"
__contributors__ = "Magdalena Musielak, Kristian Rother"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"


from unittest import main, TestCase    # for self.assertEqualItems
from moderna.analyze.ClashRecognizer import ClashRecognizer
from moderna.ModernaStructure import ModernaStructure
from Bio.PDB.PDBParser import PDBParser
from test_data import TEST_DATA_PATH

CLASH = TEST_DATA_PATH+'clash/clash.ent'
NOCLASH = TEST_DATA_PATH+'clash/no_clash.ent'
FIRST = TEST_DATA_PATH+'clash/first_model.pdb'

class FindClashesTests(TestCase):

    def setUp(self):
        self.cr = ClashRecognizer()

    def test_clash_recognition_structure(self):
        """Should find clashes in Bio.PDB.Structure object"""
        s = PDBParser().get_structure('test',CLASH)
        self.assertEqual(\
        str(self.cr.find_clashes_in_structure(s)), \
        "\
[(<Residue   G het=  resseq=3 icode= >, <Residue   G het=  resseq=4 icode= >)]\
")
        s = PDBParser().get_structure('test', NOCLASH)
        self.assertFalse(self.cr.find_clashes_in_structure(s))
        
    def test_clash_recognition_moderna_structure(self):
        """Should find clashes in ModernaStructure object"""
        s = ModernaStructure('file', CLASH)
        self.assertEqual(str(self.cr.find_clashes_in_residues(s)), \
        "[(<Residue 3 G>, <Residue 4 G>)]")
        
        s = ModernaStructure('file', NOCLASH)
        self.assertFalse(self.cr.find_clashes_in_residues(s))
        
    def test_clash_recognition_empty(self):
        """Shouldn't find any clashes in an empty list"""
        self.assertFalse(self.cr.find_clashes_in_residues([]))
        
    def test_clash_recognition_too_close_P_and_O3(self):
        """Shouldn't find any clashes if P & O3 from neighbors in 1.4A"""
        # This feature was requested by MM and KR
        s = ModernaStructure('file', FIRST)
        self.assertEqual(self.cr.find_clashes_in_residues(s), [])
        
    def test_clash_recognition(self):
        """The two example files should be distinguished."""
        self.assertEqual(\
        str(self.cr.find_clashes_in_pdb(CLASH)), \
        "[(<Residue   G het=  resseq=3 icode= >, \
<Residue   G het=  resseq=4 icode= >)]")

        self.assertFalse(\
        self.cr.find_clashes_in_pdb(NOCLASH))

    def test_two_structures(self):
        """Residue lists from two structures should be valid input."""
        # 
        s1 = PDBParser().get_structure('test', NOCLASH)
        s2 = PDBParser().get_structure('test', NOCLASH)
        c1 = s1[0].get_list()[0]
        residues1 = [r for r in c1.get_list()]
        c2 = s2[0].get_list()[0]
        residues2 = [r for r in c2.get_list()]
        # two sets of residues in the same place.
        overlapping = residues1 + residues2
        reference_results = [\
            ["<Residue   G het=  resseq=1 icode= >", \
            "<Residue   G het=  resseq=1 icode= >"], \
            ["<Residue   G het=  resseq=3 icode= >", \
            "<Residue   G het=  resseq=3 icode= >"], \
            ["<Residue   G het=  resseq=3 icode= >", \
            "<Residue   G het=  resseq=3 icode= >"], \
            ["<Residue   C het=  resseq=2 icode= >", \
            "<Residue   C het=  resseq=2 icode= >"], \
            ["<Residue   G het=  resseq=4 icode= >", \
            "<Residue   G het=  resseq=4 icode= >"], \
            ["<Residue   C het=  resseq=2 icode= >", \
            "<Residue   C het=  resseq=2 icode= >"], \
            ["<Residue   G het=  resseq=4 icode= >", \
            "<Residue   G het=  resseq=4 icode= >"], \
            ["<Residue   A het=  resseq=5 icode= >", \
            "<Residue   A het=  resseq=5 icode= >"], \
            ["<Residue   A het=  resseq=5 icode= >", \
            "<Residue   A het=  resseq=5 icode= >"], \
            ["<Residue   G het=  resseq=1 icode= >", \
            "<Residue   G het=  resseq=1 icode= >"], \
            ["<Residue   U het=  resseq=6 icode= >", \
            "<Residue   U het=  resseq=6 icode= >"], \
            ["<Residue   U het=  resseq=6 icode= >", \
            "<Residue   U het=  resseq=6 icode= >"]]

        tmp = self.cr.find_clashes_in_residues(overlapping)
        real_results = map(lambda smth: [str(smth[0]), str(smth[1])], tmp)

        # KR: replaced this for runnability without cogent
        # self.assertEqualItems(real_results, reference_results)
        for x in real_results:
            self.assert_(x in reference_results)
        
        # residues from different structures but not overlapping
        resis = residues1[:2] + residues2[-2:]
        self.assertFalse(self.cr.find_clashes_in_residues(resis))

    #
    # assumed covalent radii:
    #   c: 0.72
    #   n: 0.69
    #
    def test_cc_clash_distance(self):
        """Tests whether you can apply different radii."""
        self.assertEqual(\
        str(self.cr.find_clashes_in_pdb(TEST_DATA_PATH+'clash/cc_142.ent')),\
        "[(<Residue   G het=  resseq=1 icode= >, \
<Residue   G het=  resseq=2 icode= >)]")

        self.assertFalse\
        (self.cr.find_clashes_in_pdb(TEST_DATA_PATH+'clash/cc_145.ent'))
        
        self.assertFalse\
        (self.cr.find_clashes_in_pdb(TEST_DATA_PATH+'clash/cc_400.ent'))
        
        self.assertEqual(\
        str(self.cr.find_clashes_in_pdb(TEST_DATA_PATH+'clash/nn_137.ent')),\
        "[(<Residue   G het=  resseq=1 icode= >, \
<Residue   G het=  resseq=2 icode= >)]")

        self.assertFalse(\
        self.cr.find_clashes_in_pdb(TEST_DATA_PATH+'clash/nn_140.ent'))

        self.assertEqual(\
        str(self.cr.find_clashes_in_pdb(TEST_DATA_PATH+'clash/nc_140.ent')),\
        "[(<Residue   G het=  resseq=1 icode= >, \
<Residue   G het=  resseq=2 icode= >)]")

        self.assertFalse(\
        self.cr.find_clashes_in_pdb(TEST_DATA_PATH+'clash/nc_142.ent'))

class MemoryLeakTests(TestCase):
    
    def test_multi(self):
        "Memory should not run full - this test should not terminate!"
        print 'THIS TEST RUNS FOREVER IF OK! - WATCH MEMORY'
        s = ModernaStructure('file', FIRST)        
        while 1:
            cr = ClashRecognizer()
            cr.find_clashes_in_residues(s)
    
if __name__ == '__main__':
    main()
    
