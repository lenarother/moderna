#!/usr/bin/env python
#
# test_alignment.py
#
# unit tests for alignment
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
from moderna.sequence.RNAAlignment import RNAAlignment
from moderna.sequence.RNAAlignment import read_alignment
from moderna.sequence.ModernaSequence import Sequence
from moderna.ModelingRecipe import RecipeMaker
from moderna.Errors import AlignmentError,  AlphabetError, SequenceError
from tests.test_data import *
import os

OUT_NAME = 'out.fasta'


class ReadAlignmentTests(TestCase):
    """
    """
    def setUp(self):
        """Loads the fasta file with alignmnet."""
        self.a = read_alignment(MINI_ALIGNMENT)
        
    def tearDown(self):
        self.a = None
        if os.access(OUT_NAME, os.F_OK):
            os.remove(OUT_NAME)
        
    def test_init_string(self):
        """Alignment should be initializable by FASTA string as well."""
        a = read_alignment(FASTA_ALIGN)
        self.assertEqual(a.aligned_template_seq, Sequence('AG-CU'))
        self.assertEqual(a.aligned_target_seq, Sequence('AGCCG'))
        
    def test_init_error(self):
        """An exception should be raised if the file does not exist."""
        self.assertRaises(AlignmentError, read_alignment, 'nonexisting_file.fasta')
        
    def test_str(self):
        """Convert to string"""
        self.assertEqual(str(self.a),"""
ACUGUGAYUA[UACCU#PG
GCGGA----UUUALCUCAG
""")
        
    def test_unequal_length(self):
        """Alignment should not be initialized if sequences differ in length."""
        self.assertRaises(AlignmentError, read_alignment, LENGTH_ERROR_ALIGN)

    def test_underscore(self):
        """Alignment should not be initialized if underscore matches other than gap."""
        self.assertRaises(AlignmentError, read_alignment, UNDERSCORE_ERROR_ALIGN)


class RecipeMakerTests(TestCase):
    def setUp(self):
        """Loads the fasta file with alignmnet."""
        self.a = read_alignment(MINI_ALIGNMENT)
        self.recipe = RecipeMaker(self.a).recipe
        
    def tearDown(self):
        self.a = None
        if os.access(OUT_NAME, os.F_OK):
            os.remove(OUT_NAME)

    def test_boxes_basic(self):
        """Makes sure that the Alignment class deals in proper way with boxes"""  
        self.assertEqual(len(self.recipe.copy), 7)
        self.assertEqual(len(self.recipe.copy_backbone), 0)
        self.assertEqual(len(self.recipe.exchange), 4)
        self.assertEqual(len(self.recipe.add_modifications), 3)
        self.assertEqual(len(self.recipe.remove_modifications), 1)
        self.assertEqual(len(self.recipe.add_fragment), 1)
        self.assertEqual(len(self.recipe.add_fragment[0]), 8)
        self.assertEqual(len(self.recipe.difficult), 0)
        
    def test_boxes_nonoverlapping(self):
        """AlignmentPositionObjects should each be in only one bolignment should not be initialized if sequences differx."""
        seen = []
        boxes = [self.recipe.copy,self.recipe.copy_backbone,self.recipe.exchange,self.recipe.add_modifications,
            self.recipe.add_fragment,self.recipe.remove_modifications,self.recipe.difficult]
        for box in boxes:
            for ap in box:
                self.assertTrue(ap not in seen)
                seen.append(ap)

    def test_add_fragment_box(self):
        """Checks the list that contains positions of loops to add."""
        frags = self.recipe.add_fragment
        self.assertEqual(len(frags),1)
        fr = frags[0]
        self.assertEqual(len(fr),8)
        template_seq = ''.join([ap.template_letter and ap.template_letter.short_abbrev or '-' for ap in fr])
        target_seq = ''.join([str(ap.target_letter.short_abbrev) for ap in fr])
        self.assertEqual(template_seq,'GA----UU')
        self.assertEqual(target_seq,'GUGAYUA[')


    def test_overhangs(self):
        """Gaps at begin and end should be used only if in the template."""
        a = read_alignment(OVERHANG_ALIGN)
        recipe = ModelingRecipe(a).recipe
        self.assertEqual(len(recipe.add_fragment_5p), 1)
        self.assertEqual(len(recipe.add_fragment_5p[0]), 3)
        self.assertEqual(len(recipe.add_fragment_3p), 1)
        self.assertEqual(len(recipe.add_fragment_3p[0]), 4)
        a = read_alignment(NO_OVERHANG_ALIGN)
        recipe = ModelingRecipe(a).recipe
        self.assertEqual(len(recipe.add_fragment_5p), 0)
        self.assertEqual(len(recipe.add_fragment_3p), 0)
        
        
    def test_close_gaps(self):
        """Gaps close to each other should also work."""
        a = read_alignment(CLOSEGAPS)
        recipe = ModelingRecipe(a).recipe
        self.assertEqual(len(recipe.difficult), 0)
        self.assertEqual(len(recipe.add_fragment), 4)
        
        
    def test_short_with_insert(self):
        a = read_alignment(WITH_INSERT)
        recipe = ModelingRecipe(a).recipe
        self.assertEqual(len(recipe.difficult), 0)
        self.assertEqual(len(recipe.add_fragment), 1)
        self.assertEqual(len(recipe.add_fragment_5p), 0)
        self.assertEqual(len(recipe.add_fragment_3p), 0)

    def find_in_box(self, recipe, ap, match):
        """Checks if the given alignment position is in the given box."""
        frags = [a for b in recipe.add_fragment for a in b]
        boxes = {
            'c':[recipe.copy], 
            'b':[recipe.copy_backbone], 
            'e':[recipe.exchange], 
            'r':[recipe.remove_modifications], 
            'a':[recipe.add_modifications], 
            'C':[recipe.copy, frags], 
            'B':[recipe.copy_backbone, frags], 
            'E':[recipe.exchange, frags], 
            'R':[recipe.remove_modifications, frags], 
            'A':[recipe.add_modifications, frags], 
            'F':[frags], 
            '5':[[a for b in recipe.add_fragment_5p for a in b]], 
            '3':[[a for b in recipe.add_fragment_3p for a in b]], 
            'X':[recipe.difficult]
        }
        if match != '.':
            found = 0
            boxes_to_match = boxes.get(match, recipe.difficult)
            for box in boxes_to_match:
                for box_ap in box:
                    if ap == box_ap: 
                        found += 1
                        break
            if found == len(boxes_to_match): 
                return True
        else:
            # dot means: should be in no box.
            for box in boxes.values():
                for box_ap in box:
                    if ap == box_ap: 
                        return False
            return True
                
    def test_boxes(self):
        for seq1, seq2, matches in BOX_SAMPLES:
            align = read_alignment("> 1\n%s\n> 2\n%s\n"%(seq1, seq2))
            for ap, m in zip(align, matches):
                self.assertTrue(self.find_in_box(align, ap, m))


class AlignmentTests(TestCase):
    """"""
    # TODO: MOdeve to RNAAlignmentTests
    def setUp(self):
        """Loads the fasta file with alignmnet."""
        self.a = read_alignment(MINI_ALIGNMENT)
        
    def tearDown(self):
        self.a = None
        if os.access(OUT_NAME, os.F_OK):
            os.remove(OUT_NAME)

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
            
    def test_get_aligned_template_seq(self):
        """Should return a string with the sequence plus gaps."""
        self.assertEqual(self.a.aligned_template_seq,Sequence('GCGGA----UUUALCUCAG'))

    def test_get_aligned_target_seq(self):
        """Should return a string with the sequence plus gaps."""
        self.assertEqual(self.a.aligned_target_seq,Sequence('ACUGUGAYUA[UACCU#PG'))

    def test_load_gapstart(self):
        """Should load an alignment that starts with a gap."""
        a = read_alignment(DNA_ALIGNMENT)
        self.assertEqual(a.template_seq,Sequence('agctgccaggcaccagtg'))

    def test_new_nomenclature_alignment(self):
        a = read_alignment(NEW_NOMENCLATURE_ALIGN)
        self.assertEqual(len(a),20)
        self.assertRaises(AlignmentError,read_alignment,NEW_NOMENCLATURE_ALIGN_ERROR)
        
    def test_create_gap_begin(self):
        """Gaps at the beginning should work."""
        a = read_alignment(GAP_BEGIN_ALIGN)

    def test_remove_excess_gaps(self):
        """Should shrink adjacent gaps together."""
        ali = read_alignment(ADJACENT_GAPS)
        self.assertEqual(str(ali),SHRUNK_GAPS)

    def test_remove_excess_gaps_underscore(self):
        """Should not shrink underscore positions."""
        ali = read_alignment(ADJACENT_GAPS_US)
        self.assertEqual(str(ali),SHRUNK_GAPS_US)

    def test_complicated_trna_alignment(self):
        ali = read_alignment("""> 1j2b_D.pdb
GGGCCCGUGGUCUAGUU_--------G-ACGCC-GCCC-UUACGAGGCGGAG-----------G-UC-CGGGGUUC-A--AG--U-C-CC--C-G-CG-GGCCCA-CCA
> 2du6_D.pdb
--GCCAGGGUGGCAGA---GGGGCUUU-GCGGC-GGAC-UCUAGAUCCGCUUU----------A--C-CCCGGUUC-G--AA--U-C-CG--G-G-CC-CUG-GC----
""")
        expected = """
GGGCCCGUGGUCUAGUU_------GACGCCGCCCUUACGAGGCGGAG-GUCCGGGGUUCAAGUCCCCGCGGGCCCACCA
--GCCAGGGUGGCAGAG-GGGCUUUGCGGCGGACUCUAGAUCCGCUUUA-CCCCGGUUCGAAUCCGGGCCCUG-GC---
"""
        self.assertEqual(str(ali), expected)



BOX_SAMPLES = [
    (
    "AGCU", 
    "AGCU", 
    "cccc"
    ), 
    (
    "AG#U", 
    "U7CU", 
    "erac"
    ), 
    (
    "AAAAAAAAA", 
    "--AAAAA--", 
    "55ccccc33"
    ), 
    (
    "--AAAAA--", 
    "AAAAAAAAA", 
    "..ccccc.."
    ), 
    (
    "AAA----AAAAA", 
    "AAAAAAAAAAAA", 
    "ccF....Fcccc"
    ), 
    (
    "AAAAAAAAAAAA", 
    "AAA----AAAAA", 
    "cCCFFFFCCccc"
    ), 
    # now the nasty cases
    (
    "AAA--------------------------AAAAA", 
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
    "ccF..........................Fcccc"
    ), 
    (
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
    "AAA--------------------------AAAAA", 
    "ccFFFFFFFFFFFFFFFFFFFFFFFFFFFFcccc"
    ), 
    (
    "AAA--------------------------AAAAAA", 
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-AAAA", 
    "ccF..........................FFFccc"
    ), 
    (
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-AAAA", 
    "AAA--------------------------AAAAAA", 
    "ccFFFFFFFFFFFFFFFFFFFFFFFFFFFF.Fccc"
    ), 
    (
    "AAAAAAAAAAAAAAAAAAAAAA", 
    "AAA-----AA------AAAAAA", 
    "ccFFFFFFFFFFFFFFFccccc"
    ), 
    (
    "AAA-----AA------AAAAAA", 
    "AAAAAAAAAAAAAAAAAAAAAA", 
    "ccF.....FF......Fccccc"
    ), 
    (
    "AAA-----A------AAAAAA", 
    "AAAAAAAAAAAAAAAAAAAAA", 
    "ccF.....F......Fccccc"
    ), 
    (
    "AAAAAAAAAAAAAAAAAAAAA", 
    "AAA-----A------AAAAAA", 
    "ccFFFFFFFFFFFFFFccccc"
    ),
    (
    "A-AAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAA",
    "5.ccccccccccccccccccc"
    ),
    (
    "AAAAAAAAAAAAAAAAAAAAA",
    "A-AAAAAAAAAAAAAAAAAAA",
    "55ccccccccccccccccccc"
    ),
    (
    "AAAAAAAAAAAAAAAAAAA-A",
    "AAAAAAAAAAAAAAAAAAAAA",
    "ccccccccccccccccccc.3"
    ),
    (
    "AAAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAA-A",
    "ccccccccccccccccccc33"
    ),
    (
    "A-----AAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAA",
    "5.....ccccccccccccccc"
    ),
    (
    "AAAAAAAAAAAAAAAAAAAAA",
    "A-----AAAAAAAAAAAAAAA",
    "555555ccccccccccccccc"
    ),
    (
    "AAAAAAAAAAAAAAA-----A",
    "AAAAAAAAAAAAAAAAAAAAA",
    "ccccccccccccccc.....3"
    ),
    (
    "AAAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAA-----A",
    "ccccccccccccccc333333"
    ),

]    

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

UNDERSCORE_ERROR_ALIGN = """>target
CCCCC
> template
CG_CG
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


if __name__ == '__main__':
    main()
  

