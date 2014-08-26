#!/usr/bin/env python
#
# test_alignment.py
#
# unit tests for base exchange functionality
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
from moderna.sequence.RNAAlignment import RNAAlignment, RNAAlignmentParser
from moderna.sequence.ModernaSequence import Sequence
from moderna.util.Errors import AlignmentError,  AlphabetError, SequenceError
from tests.test_data import *
import os

OUT_NAME = 'out.fasta'


class RNAAlignmentTests(TestCase):
    """
    """
    def setUp(self):
        """Loads the fasta file with alignmnet."""
        self.parser = RNAAlignmentParser()
        self.a = self.parser.get_alignment(MINI_ALIGNMENT)
        
    def tearDown(self):
        self.a = None
        if os.access(OUT_NAME, os.F_OK):
            os.remove(OUT_NAME)
        
    def test_init_string(self):
        """Alignment should be initializable by FASTA string as well."""
        a = self.parser.get_alignment(FASTA_ALIGN)
        self.assertEqual(a.aligned_template_seq, Sequence('AG-CU'))
        self.assertEqual(a.aligned_target_seq, Sequence('AGCCG'))
        
    def test_init_error(self):
        """An exception should be raised if the file does not exist."""
        self.assertRaises(AlignmentError, self.parser.get_alignment_from_file, 'nonexisting_file.fasta')
        
    def test_str(self):
        """Convert to string"""
        self.assertEqual(str(self.a),"""
ACUGUGAYUA[UACCU#PG
GCGGA----UUUALCUCAG
""")
        
    def test_unequal_length(self):
        """Alignment should not be initialized if sequences differ in length."""
        self.assertRaises(AlignmentError, self.parser.get_alignment, LENGTH_ERROR_ALIGN)
        
    def test_getitem(self):
        """Alignment positions should be indexed starting with 1."""
        self.assertEqual(self.a[1].template_letter.short_abbrev,'G')
        self.assertEqual(self.a[6].template_letter,None)
        self.assertEqual(self.a[10].template_letter.short_abbrev,'U')
        self.assertEqual(self.a[14].template_letter.short_abbrev,'L')
        self.assertEqual(self.a[1].target_letter.short_abbrev,'A')
        self.assertEqual(self.a[6].target_letter.short_abbrev,'G')
        self.assertEqual(self.a[10].target_letter.short_abbrev,'A')
        self.assertEqual(self.a[14].target_letter.short_abbrev,'C')

    def test_position_indices(self):
        """Alignment position indices need to be consistent."""
        for i in range(1,10):
            self.assertEqual(i,self.a[i].alignment_position)        

    def test_set_aligned_sequences(self):
        """Both sequences can be set"""
        self.a.set_aligned_sequences(Sequence('AAUU'), Sequence('GGCC'))
        self.assertEqual(self.a.aligned_target_seq, Sequence('AAUU'))
        self.assertEqual(self.a.aligned_template_seq, Sequence('GGCC'))

    def test_set_aligned_sequences_fail(self):
        """Both sequences need to have equal length"""
        self.assertRaises(AlignmentError, self.a.set_aligned_sequences, Sequence('AAUU'), Sequence('GGC'))

    def test_get_aligned_template_seq(self):
        """Should return a string with the sequence plus gaps."""
        self.assertEqual(self.a.aligned_template_seq, Sequence('GCGGA----UUUALCUCAG'))

    def test_get_aligned_target_seq(self):
        """Should return a string with the sequence plus gaps."""
        self.assertEqual(self.a.aligned_target_seq, Sequence('ACUGUGAYUA[UACCU#PG'))

    def test_load_gapstart(self):
        """Should load an alignment that starts with a gap."""
        a = self.parser.get_alignment(DNA_ALIGNMENT)
        self.assertEqual(a.template_seq,Sequence('agctgccaggcaccagtg'))

    def test_new_nomenclature_alignment(self):
        a = self.parser.get_alignment(NEW_NOMENCLATURE_ALIGN)
        self.assertEqual(len(a),20)
        self.assertRaises(AlignmentError,self.parser.get_alignment,NEW_NOMENCLATURE_ALIGN_ERROR)
        
    def test_create_gap_begin(self):
        """Gaps at the beginning should work."""
        a = self.parser.get_alignment(GAP_BEGIN_ALIGN)
        
    def test_remove_excess_gaps(self):
        """Should shrink adjacent gaps together."""
        ali = self.parser.get_alignment(ADJACENT_GAPS)
        self.assertEqual(str(ali),SHRUNK_GAPS)

    def test_remove_excess_gaps_underscore(self):
        """Should not shrink underscore positions."""
        ali = self.parser.get_alignment(ADJACENT_GAPS_US)
        self.assertEqual(str(ali),SHRUNK_GAPS_US)

    def test_complicated_trna_alignment(self):
        ali = self.parser.get_alignment(COMPLICATED_ALIGN)
        expected = """
GGGCCCGUGGUCUAGUU_------GACGCCGCCCUUACGAGGCGGAG-GUCCGGGGUUCAAGUCCCCGCGGGCCCACCA
--GCCAGGGUGGCAGAG-GGGCUUUGCGGCGGACUCUAGAUCCGCUUUA-CCCCGGUUCGAAUCCGGGCCCUG-GC---
"""
        self.assertEqual(str(ali), expected)

    def test_write_fasta(self):
        """Writes FASTA files"""
        a = self.parser.get_alignment(NEW_NOMENCLATURE_ALIGN)
        a.write_fasta_file(OUT_NAME)
        self.assertTrue(OUT_NAME, os.F_OK)
        a2 = self.parser.get_alignment_from_file(OUT_NAME)
        self.assertEqual(a2.target_seq, a.target_seq)
        self.assertEqual(a2.template_seq, a.template_seq)
        
    def test_write_fasta_nomenclature(self):
        """Writes FASTA files with new nomenclature"""
        a = self.parser.get_alignment(NEW_NOMENCLATURE_ALIGN)
        a.write_fasta_file(OUT_NAME, new_nomenclature=True)
        self.assertTrue(OUT_NAME, os.F_OK)
        a2 = self.parser.get_alignment_from_file(OUT_NAME)
        self.assertEqual(a2.target_seq.seq_new_notation, a.target_seq.seq_new_notation)
        self.assertEqual(a2.template_seq.seq_new_notation, a.template_seq.seq_new_notation)
        
    def test_calculate_dissimilarity(self):
        """Should return dissimilarity between target and template."""
        self.assertAlmostEqual(0.60526315789473684,  self.a.calculate_dissimilarity(), 5)
        ali = self.parser.get_alignment(ALIGN_1B23_1QF6)
        self.assertAlmostEqual(0.41973684210526313,  ali.calculate_dissimilarity(), 5)
        ali = self.parser.get_alignment(ALIGN_TARGET_GAP)
        self.assertAlmostEqual(0.76666666666666672,  ali.calculate_dissimilarity(), 5)
        self.assertAlmostEqual(ali.calculate_dissimilarity(),  (1.0-ali.calculate_similarity()), 5)
        ali = self.parser.get_alignment(DNA_ALIGNMENT)
        self.assertAlmostEqual(0.80952380952380953,  ali.calculate_dissimilarity(), 5)
        ali = self.parser.get_alignment(SIMILARITY_ALIGNMENT)
        self.assertAlmostEqual(0.50769230769230766,  ali.calculate_dissimilarity(), 5)

    def test_calculate_similarity(self):
        """Should return similarity between target and template."""
        self.assertAlmostEqual(0.39473684210526316, self.a.calculate_similarity(), 5)
        ali = self.parser.get_alignment(ALIGN_1B23_1QF6)
        self.assertAlmostEqual(0.58026315789473681,  ali.calculate_similarity(), 5)
        ali = self.parser.get_alignment(ALIGN_TARGET_GAP)
        self.assertAlmostEqual(0.23333333333333334,  ali.calculate_similarity(), 5)
        ali = self.parser.get_alignment(DNA_ALIGNMENT)
        self.assertAlmostEqual(0.19047619047619047,  ali.calculate_similarity(), 5)
        ali = self.parser.get_alignment(SIMILARITY_ALIGNMENT)
        self.assertAlmostEqual(0.49230769230769228,  ali.calculate_similarity(), 5)
        
    def test_calculate_identity(self):
        """Should return identity between target and template."""
        self.assertAlmostEqual(0.36842105263157893,  self.a.calculate_identity(), 5)
        ali = self.parser.get_alignment(ALIGN_1B23_1QF6)
        self.assertAlmostEqual(0.48684210526315791,  ali.calculate_identity(), 5)
        ali = self.parser.get_alignment(ALIGN_TARGET_GAP)
        self.assertAlmostEqual(0.066666666666666666,  ali.calculate_identity(), 5)
        ali = self.parser.get_alignment(DNA_ALIGNMENT)
        self.assertAlmostEqual(0.11904761904761904,  ali.calculate_identity(), 5)
        ali = self.parser.get_alignment(SIMILARITY_ALIGNMENT)
        self.assertAlmostEqual(0.20512820512820512,  ali.calculate_identity(), 5)
        
    

class RNAAlignmentParserTests(TestCase):
    
    def setUp(self):
        self.parser = RNAAlignmentParser()
        
    def test_get_alignment(self):
        result = self.parser.get_alignment(MINI_ALIGNMENT)
        self.assertEqual(type(result), RNAAlignment)
    


SHRUNK_GAPS = """
-GGUUUCCCAAA
AAAUUUCCCGGG
"""

SHRUNK_GAPS_US = """
-GGUUU-CCCAAA_--CC
AAAUUU_CCCGGG-CCCC
"""


ADJACENT_GAPS = """> second
---GGU--UUCCCAAA---
> first
AAA--UUU--CCC---GGG
"""

ADJACENT_GAPS_US = """> second
---GGU---UUCCCAAA_-----CC
> first
AAA--UUU_--CCC----GGGCCCC
"""

CLOSEGAPS = '''> target
GGGAUAGUUCCAGACCCCBCCCCU#A
> template
GGGA-AG--CCAGA----B----U#A
'''

WITH_INSERT = '''> target
GG-AAAACC
> template
GG_----CC
'''


FASTA_ALIGN = """> target
AGCCG
> template
AG-CU
"""

LENGTH_ERROR_ALIGN = """> target
AGCCG
> template
AG-CUU
"""

GAP_BEGIN_ALIGN = """> 2
-AAAA-----------AAAAAAAAAAA
> 1
A-AAAAAAAAAAAAAA-----------
"""

OVERHANG_ALIGN = """> 2
AAAAAAAAAAAAAAAAAA
> 1
---AAAAAAAAAAA----
"""

NO_OVERHANG_ALIGN = """> 2
---AAAAAAAAAAA----
> 1
AAAAAAAAAAAAAAAAAA
"""

NEW_NOMENCLATURE_ALIGN = """> model sequence for testing
ACUGUGAYUA[U   UACCU#PG
> uses number-encoding of modifications
GCGGA----UUU015UALCUCAG
"""

NEW_NOMENCLATURE_ALIGN_ERROR = """> model sequence for testing
ACUGUGAYUA[U--UUACCU#PG
> uses number-encoding of modifications
GCGGA----UUU015UALCUCAG
"""

COMPLICATED_ALIGN = """> 1j2b_D.pdb
GGGCCCGUGGUCUAGUU_--------G-ACGCC-GCCC-UUACGAGGCGGAG-----------G-UC-CGGGGUUC-A--AG--U-C-CC--C-G-CG-GGCCCA-CCA
> 2du6_D.pdb
--GCCAGGGUGGCAGA---GGGGCUUU-GCGGC-GGAC-UCUAGAUCCGCUUU----------A--C-CCCGGUUC-G--AA--U-C-CG--G-G-CC-CUG-GC----
"""

if __name__ == '__main__':
    main()
  

