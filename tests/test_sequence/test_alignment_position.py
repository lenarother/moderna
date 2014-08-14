#!/usr/bin/env python
#
# test_alignment_position.py
#
# unit tests for SequenceAlignment.AlignmentPosition
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
from moderna.RNAAlignment import AlignmentPosition
from moderna.Constants import ANY_RESIDUE
from moderna.ModernaAlphabet import AlphabetEntry
from test_data import *


class AlignmentPositionTests(TestCase):

    def setUp(self):
        """Creates some alphabet entries"""

    def test_attributes(self):
        """Init should set all four attributes."""
        ap = AlignmentPosition(TAR_P,'tarL',TEM_P,'temL')
        self.assertEqual(ap.target_position,TAR_P)
        self.assertEqual(ap.template_position,TEM_P)
        self.assertEqual(ap.target_letter,'tarL')
        self.assertEqual(ap.template_letter,'temL')

    def test_is_identical(self):
        """identical letters and wildcard combinations should return True"""
        # None is a gap
        expected = [True, False, False, False, False, False, False, False, True]
        for data, result in zip(CASES, expected):
            tar_let, tem_let = data
            ap = AlignmentPosition(TAR_P,tar_let,TEM_P,tem_let)
            self.assertEqual(ap.is_identical(),result)

    def test_is_different(self):
        """Checks which AP's are different."""
        expected = [False, True, True, True, True, True, True, True, False]
        # loop all test cases
        for data, result in zip(CASES, expected):
            tar_let, tem_let = data
            ap = AlignmentPosition(TAR_P,tar_let,TEM_P,tem_let)
            self.assertEqual(ap.is_different(),result)
            
    def test_is_gap(self):
        # gap mode
        expected = [False, False, True, True, False, False, True, True, True]
        # loop all test cases
        for data, result in zip(CASES, expected):
            tar_let, tem_let = data
            ap = AlignmentPosition(TAR_P,tar_let,TEM_P,tem_let)
            self.assertEqual(ap.has_gap(),result)
            
    def test_is_gap_in_template(self):
        # gap in template mode
        expected = [False, False, True, False, False, False, False, True, True]
        # loop all test cases
        for data, result in zip(CASES, expected):
            tar_let, tem_let = data
            ap = AlignmentPosition(TAR_P,tar_let,TEM_P,tem_let)
            self.assertEqual(ap.has_template_gap(),result)
            
    def test_is_gap_in_target(self):
        # gap in target mode
        expected = [False, False, False, True, False, False, True, False, True]
        # loop all test cases
        for data, result in zip(CASES, expected):
            tar_let, tem_let = data
            ap = AlignmentPosition(TAR_P,tar_let,TEM_P,tem_let)
            self.assertEqual(ap.has_target_gap(),result)
    
    def test_is_mismatch(self):
        # mismatch mode
        expected = [False, True, False, False, False, False, False, False, False]
        # loop all test cases
        for data, result in zip(CASES, expected):
            tar_let, tem_let = data
            ap = AlignmentPosition(TAR_P,tar_let,TEM_P,tem_let)
            self.assertEqual(ap.is_mismatch(),result)
    
    def test_is_unidentified(self):
        # is unidentified mode
        expected = [False, False, False, False, True, True, False, False, False]
        # loop all test cases
        for data, result in zip(CASES, expected):
            tar_let, tem_let = data
            ap = AlignmentPosition(TAR_P,tar_let,TEM_P,tem_let)
            self.assertEqual(ap.is_unidentified(),result)

# example AlphabetEntries
AE_M = AlphabetEntry('mmm','M','MAM','magdalena','M','10M')
AE_K = AlphabetEntry('kkk','K','KAK','kristian','K','11K')
AE_R = AlphabetEntry('rrr','R','RAR','rother','R','12R')
AE_ANY = AlphabetEntry('X','.','X','Unknown','X','UNK')

# example target and template positions
TAR_P = 10
TEM_P = 20

# example data for AlignmentPositions
CASES = [
        (AE_M, AE_M), 
        (AE_K, AE_R), 
        (AE_M, None), 
        (None, AE_R), 
        (AE_M, AE_ANY), 
        (AE_ANY, AE_ANY), 
        (None, AE_ANY), 
        (AE_ANY, None), 
        (None, None) # stupid case but should be defined.
        ]

if __name__ == '__main__':
    main()
