#!/usr/bin/env python
#
# test_sequence.py
#
# unit tests for sequence class
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
from moderna.ModernaSequence import Sequence
from moderna.ModernaAlphabet import AlphabetEntry
from moderna.Errors import SequenceError, AlphabetError

class SequenceTests(TestCase):
    """
    """
    def setUp(self):
        """Defines some simple cases."""
        self.normal = 'AGCU'
        self.longer = 'AAAAAAAAAAA'
        self.modified = 'CCCCDDDGGG'
        self.unmodified = 'CCCCUUUGGG'        
        self.different = 'CCCC777GGG'

    def test_init(self):
        """Attributes of seq should work."""
        s = Sequence(self.normal)
        self.assertEqual(s.seq_with_modifications,'AGCU')
        self.assertEqual(s.seq_without_modifications,'AGCU')
        s = Sequence(self.modified)
        self.assertEqual(s.seq_with_modifications,'CCCCDDDGGG')
        self.assertEqual(s.seq_without_modifications,'CCCCUUUGGG')

    def test_create_sequence(self):
        """All base characters should be allowed in a sequence."""
        s = "ACGUNRYagct"
        seq = Sequence(s)
        self.assertEqual(seq.seq_with_modifications,s)
        
    def test_create_sequence_mods(self):
        """All modification characters should be allowed in a sequence."""
        s = "P[]*#37T<>"
        # KR: add remaining modifications here
        seq = Sequence(s)
        self.assertEqual(seq.seq_with_modifications,s)
        
    def test_create_sequence_error(self):
        """Non-base characters should raise an error.""" 
        for char in " ":
            self.assertRaises(SequenceError,Sequence,"AC"+char+"GU")
        for char in "n":
            self.assertRaises(AlphabetError,Sequence,"AC"+char+"GU")
        
    def test_equal_unmodified(self):
        """unmodified sequences should be comparable."""
        s1 = Sequence(self.normal)
        s2 = Sequence(self.normal)
        s3 = Sequence(self.longer)
        self.assertEqual(s1,s2)
        self.assertNotEqual(s1,s3)

    def test_equal_modified(self):
        """modified sequences should NOT match."""
        s1 = Sequence(self.modified)
        s2 = Sequence(self.modified)
        s3 = Sequence(self.unmodified)
        s4 = Sequence(self.different)
        self.assertEqual(s1,s2)
        self.assertNotEqual(s1,s3)
        self.assertNotEqual(s1,s4)

    def test_similar(self):
        """modified sequences are similar to their unmodified counterparts."""
        s1 = Sequence(self.modified)
        s2 = Sequence(self.modified)
        s3 = Sequence(self.unmodified)
        s4 = Sequence(self.different)
        self.assertTrue(s1.similar_to(s2))
        self.assertTrue(s1.similar_to(s3))
        self.assertTrue(s3.similar_to(s1))
        self.assertTrue(s3.similar_to(s3))
        self.assertFalse(s3.similar_to(s4))
        self.assertFalse(s1.similar_to(s4))

    def test_standard_chars(self):
        """All normal characters should work."""
        chars = "ACGUacgt"
        s = Sequence(chars)
        self.assertEqual(s.seq_with_modifications, chars)

    def test_all_modifications(self):
        """All modification characters should work."""
        all_mods = """/"+*=6E[:IO^`%BM?'}>KL#R|7(Q89YW{2J4&1S3V5!$X,)~DP]ZTF\\"""
        s = Sequence(all_mods)
        self.assertEqual(s.seq_with_modifications, all_mods)

    def test_new_modifications(self):
        """Modifications using the new nomenclature should work."""
        mods = "AAA010AA017UAA"
        s = Sequence(mods)
        self.assertEqual(len(s), 8)
        self.assertEqual(s.seq_without_modifications, "AAAAAUAA")
        self.assertRaises(SequenceError,Sequence,'AAA0AAAAAAA')
        self.assertRaises(SequenceError,Sequence,'AAA09A1AAA')
        self.assertRaises(AlphabetError,Sequence,'AAA099A1AAA')

    def test_new_nomenclature_attributes(self):
        """Should represent new nomenclature in sequence"""
        s = Sequence("GCGGAUUU015UALCUCAG")
        self.assertEqual(s.seq_without_modifications, "GCGGAUUUUAGCUCAG")
        self.assertEqual(s.seq_with_modifications, "GCGGAUUUxALCUCAG")
        self.assertEqual(s.seq_new_notation, "GCGGAUUU015UA002GCUCAG")

    def test_exclude_space(self):
        """Spaces in sequences are disallowed."""
        self.assertRaises(SequenceError,Sequence,'AAA 001AAA')
        self.assertRaises(SequenceError,Sequence,' AAA  AAAAAA')

    def test_special_chars(self):
        """Special characters like gap symbols should work."""
        # gap and unknown symbols
        s = Sequence("._-")
        self.assertEqual(s.seq_with_modifications, "._-")
        # wildcard bases
        s = Sequence("H;<N")
        self.assertEqual(s.seq_with_modifications, "H;<N")

    def test_iterate_alphabet_list(self):
        s = Sequence("AGCU001A")
        alph_entries = [ae for ae in s.seq_alphabet_list]
        self.assertEqual(len(alph_entries),5)
        self.assertTrue(isinstance(alph_entries[0],AlphabetEntry)) # ?
        
    def test_iterate(self):
        """Sequences should be iterable"""
        s = Sequence("ADADA")
        positions = [ae for ae in s]
        self.assertEqual(len(positions),5)
        self.assertTrue(isinstance(positions[0],AlphabetEntry)) 
        
    def test_iterate_twice(self):
        """Sequences should be iterable more than once"""
        s = Sequence("ADADA")
        positions = [ae for ae in s]
        positions2 = [ae for ae in s]
        self.assertEqual(len(positions2),5)
        self.assertTrue(isinstance(positions2[0],AlphabetEntry)) 


if __name__ == '__main__':
        main()
  

