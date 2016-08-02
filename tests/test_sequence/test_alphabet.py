#!/usr/bin/env python
"""
unit tests for Alphabet
"""

from unittest import main, TestCase
from moderna.sequence.ModernaAlphabet import Alphabet, AlphabetEntry
from moderna.Constants import UNKNOWN_RESIDUE_SHORT, \
    UNKNOWN_RESIDUE_ONELETTER,  RESIDUE_WITHOUT_ONE_LETTER_ABBREV, \
    MODIFICATION_NAMES_TABLE_PATH
from tests.test_data import *


class AlphabetTests(TestCase):
    """
    Tests the functionality of Alphabet as in the class
    description.
    KR won't write any further tests for single methods to
    allow the API to be flexible.
    """
    def test_init(self):
        """init should result in a dictionary with the test file contents."""
        a = Alphabet(TEST_MODIFICATION_TABLE)
        self.assertEqual(len(a.keys()), 2)
        self.assertTrue('first' in a)
        self.assertTrue('second' in a)
        self.assertEqual(a['first'].pdb_abbrev, 'UNK')
        self.assertEqual(a['first'].original_base, 'A')
        self.assertEqual(a['first'].long_abbrev, 'first')
        self.assertEqual(a['first'].short_abbrev, 'F')
        self.assertEqual(a['first'].full_name, 'first_full_name')
        self.assertEqual(a['second'].pdb_abbrev, 'UNK')
        self.assertEqual(a['second'].original_base, 'B')
        self.assertEqual(a['second'].long_abbrev, 'second')
        self.assertEqual(a['second'].short_abbrev, 'S')
        self.assertEqual(a['second'].full_name, 'second_full_name')

    def test_second_dict(self):
        """The reverse modification_abbrev->original_base dict should work."""
        a = Alphabet(TEST_MODIFICATION_TABLE)
        self.assertEqual(a._short_original['F'].original_base, 'A')
        self.assertEqual(a._short_original['S'].original_base, 'B')
        self.assertRaises(KeyError,a._short_original.__getitem__, 'A')

    def test_unknown_residue(self):
        """The modification_abbrev dict should also work for unspecified residues."""
        a = Alphabet(MODIFICATION_NAMES_TABLE_PATH)
        self.assertEqual(a.get_short_original(UNKNOWN_RESIDUE_ONELETTER).long_abbrev,UNKNOWN_RESIDUE_SHORT)
    
    def test_get_short_original_special(self):
        a = Alphabet(MODIFICATION_NAMES_TABLE_PATH)
        self.assertEqual(a.get_short_original(RESIDUE_WITHOUT_ONE_LETTER_ABBREV),a[RESIDUE_WITHOUT_ONE_LETTER_ABBREV])
        
    def test_strange_symbols(self):
        """The modification_abbrev dict should also work for misc symbols."""
        a = Alphabet(MODIFICATION_NAMES_TABLE_PATH)
        for symbol in "_-.@":
            self.assertTrue(a.get_short_original(symbol))
        
    #def test_html_signs_disabled(self):
    #    """Requesting the degenerate HTML replacement symbol as one-letter abbrev should not work."""
    #    a = Alphabet(MODIFICATION_TABLE)
    #    self.assertRaises(AlphabetError,a.get_short_original,RESIDUE_WITHOUT_ONE_LETTER_ABBREV)


class AlphabetEntryTests(TestCase):
    """
    Tests the functionality of AlphabetEntry as in the class
    description.
    KR won't write any further tests for single methods to
    allow the API to be flexible.
    """
    def test_equal(self):
        """aaa"""
        ae1=AlphabetEntry('m5C','?','UNK','5-methylcytidine','C','5C')
        ae2=AlphabetEntry('m5C','?','UNK','5-methylcytidine','C','5C')
        ae3=AlphabetEntry('ac4C','M','4AC','N4-acetylcytidine','C','6C')
        self.assertEqual(ae1,ae2)
        self.assertEqual(ae3,ae3)
        self.assertNotEqual(ae1,ae3)
        self.assertNotEqual(ae2,ae3)
        
    def test_true(self):
        """AlphabetEntries are always True"""
        ae=AlphabetEntry('m5C','?','UNK','5-methylcytidine','C','5C')
        self.assertTrue(ae)


if __name__ == '__main__':
    main()
    
