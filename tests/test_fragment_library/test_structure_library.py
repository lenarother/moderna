#!/usr/bin/env python
#
# test_moderna_structure.py
#
# unit tests for StructureLibrary class
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

from unittest import main, TestCase
from moderna.ModernaStructure import ModernaStructure, ModernaResidue
from moderna.sequence.ModernaSequence import Sequence
from moderna.fragment_library.StructureLibrary import StructureLibrary
from moderna.Constants import LIR_DIRECTORY_PATH
from test_data import *

class StructureLibraryTests(TestCase):

    def test_init(self):
        """Makes sure there is no error during initialization."""
        sl = StructureLibrary(RNA_PATH)

    def test_get_structure(self):
        """Generates ModernaStructure objects"""
        sl = StructureLibrary(RNA_PATH)
        s = sl.get_structure('minirna_14res_ss.pdb', MINI_TEMPLATE_CHAIN_NAME)
        self.assertTrue(isinstance(s, ModernaStructure))
        self.assertEqual(s.get_sequence(),Sequence("GCGGAUUUALCUCAG"))

    def test_get_structure_multiple(self):
        """Generates the same ModernaStructure object repeatedly"""
        sl = StructureLibrary(RNA_PATH)
        for i in range(100):
            s = sl.get_structure('minirna_14res_ss.pdb', MINI_TEMPLATE_CHAIN_NAME)
            self.assertEqual(s.get_sequence(),Sequence("GCGGAUUUALCUCAG"))
            
    def test_get_multiple_manipulated(self):
        """The retrieved instances may be manipulated"""
        sl = StructureLibrary(RNA_PATH)
        for i in range(10):
            s = sl.get_structure('minirna_14res_ss.pdb', MINI_TEMPLATE_CHAIN_NAME)
            self.assertEqual(s.get_sequence(),Sequence("GCGGAUUUALCUCAG"))
            s.remove_residue(str(i+1))
    
    def test_find_resi_in_lines(self):
        """Should return correct indices"""
        sl = StructureLibrary(RNA_PATH)
        lines = open(RNA_ATOMSONLY).readlines()
        begin, end = sl.find_resi_in_lines(lines, '4')
        self.assertEqual(begin, 67)
        self.assertEqual(end, 90)
        
    def test_get_structure_part(self):
        """should read small pieces"""
        sl = StructureLibrary(RNA_PATH)
        struc = sl.get_structure_part('minirna_14res_ss.pdb', 'A', '3', '7')
        self.assertTrue(isinstance(struc, ModernaStructure))
        self.assertEqual(len(struc), 5)
        
    def test_read_fragment_from_ribosome(self):
        """Should read pieces of ribosomes quickly"""
        for i in range(10):
            sl = StructureLibrary(LIR_DIRECTORY_PATH)
            ribo = 'rr0082H09.pdb'
            struc = sl.get_structure_part(ribo, '0', '107', '109')
            self.assertTrue(isinstance(struc, ModernaStructure))
            self.assertEqual(len(struc), 3)
            resids = [r.identifier for r in struc]
            self.assertEqual(resids, ['107', '108', '109'])

if __name__ == '__main__':
    import cProfile
    #cProfile.run("main()")
    main()
    
