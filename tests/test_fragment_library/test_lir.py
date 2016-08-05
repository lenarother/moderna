#!/usr/bin/env python
"""
Unit tests for linker database entries
"""

from unittest import main, TestCase
from moderna.fragment_library.LIR import LirRecord,Lir
from moderna.ModernaStructure import ModernaStructure
from moderna.sequence.ModernaSequence import Sequence
from moderna.tests.test_data import *
import os, math

class LirRecordTests(TestCase):

    def setUp(self):
        self.struc = ModernaStructure('file',MINI_TEMPLATE)
        self.lr = LirRecord( 3, '1ehz', 'A', 4, 8, Sequence('GCC'),Sequence('GGCCG'), "(...)", \
                  1.23, 4.56, 7.89, 0.12, 3.45, 6.78, 9.10, \
                  1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9)

    def test_init(self):
        """Attributes should be present."""
        self.assertEqual(self.lr.fr_length, 3)
        self.assertEqual(self.lr.structure, '1ehz')
        self.assertEqual(self.lr.chain,'A')
        self.assertEqual(self.lr.preceding_residue, 4)
        self.assertEqual(self.lr.following_residue, 8)
        self.assertEqual(self.lr.sequence, Sequence('GCC'))
        self.assertEqual(self.lr.sequence_anchor, Sequence('GGCCG'))
        self.assertEqual(self.lr.secstruc, "(...)")
        self.assertEqual(self.lr.x, 1.23)
        self.assertEqual(self.lr.y, 4.56)
        self.assertEqual(self.lr.dist_anchor, 7.89)
        self.assertEqual(self.lr.beta, 0.12)
        self.assertEqual(self.lr.gamma, 3.45)
        self.assertEqual(self.lr.omega5, 6.78)
        self.assertEqual(self.lr.omega3, 9.10)
        dd = [d for d in self.lr.distances]
        self.assertEqual(dd, [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9])

    def test_get_list(self):
        """Should return a list with all attributes."""
        l = self.lr.get_list()
        self.assertEqual(l, [3, '1ehz', 'A', 4, 8, Sequence('GCC'),Sequence('GGCCG'), "(...)", \
                              1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9])

    def test_get_txt(self):
        """Should return all attributes as text."""
        s = self.lr.get_txt('##')
        expected = "3##1ehz##A##4##8##GCC##GGCCG##(...)##1.1##2.2##3.3##4.4##5.5##6.6##7.7##8.8##9.9"
        self.assertEqual(s, expected)
            
            

if __name__ == '__main__':
    main()
  
