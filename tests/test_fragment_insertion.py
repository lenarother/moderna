#!/usr/bin/env python
#
# test_fragment_insertion.py
#
# unit tests for FragmentInserter
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
from moderna.ModernaFragment import ModernaFragment5,  \
    ModernaFragment3, ModernaFragment53
from moderna.FragmentInsertion import FragmentInserter
from moderna.RNAModel import RnaModel
from moderna.ModernaStructure import ModernaStructure
from moderna.sequence.ModernaSequence import Sequence
from moderna.Errors import ModernaFragmentError
from test_data import *

class FragmentInserterTests(TestCase):
    """
    Tests for the FragmentInserter class
    """
    def setUp(self):
        self.m = RnaModel(data_type='file',data=MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
        self.s1 = ModernaStructure('file', SMALL_FRAGMENT, seq=Sequence("GCGG"))
        #self.s2 = ModernaStructure('file',SMALL_FRAGMENT, seq=Sequence("GCGG"))
        #self.s3 = ModernaStructure('file', MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))

    def test_insert(self):
        """Inserts fragment into a model"""
        f1 = ModernaFragment53(self.s1, anchor5=self.m['10'],anchor3=self.m['13'])
        finsert = FragmentInserter()
        finsert.insert_fragment(f1, self.m)
        self.assertEqual(self.m.get_sequence(),Sequence("GCGGAUUUALCGCAG"))


if __name__ == '__main__':
    main()
