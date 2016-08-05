
from unittest import TestCase, main
from moderna.ModelingRecipe import RecipeMaker
from moderna.sequence.RNAAlignment import read_alignment
import os

from moderna.tests.test_data import MINI_ALIGNMENT

class RecipeMakerTests(TestCase):
    def setUp(self):
        """Loads the fasta file with alignmnet."""
        self.a = read_alignment(MINI_ALIGNMENT)
        self.recipe = RecipeMaker(self.a).recipe

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


    def test_overhang(self):
        """Gaps at begin and end should be used only if in the template."""
        a = read_alignment(OVERHANG_ALIGN)
        recipe = RecipeMaker(a).recipe
        self.assertEqual(len(recipe.add_fragment_5p), 1)
        self.assertEqual(len(recipe.add_fragment_5p[0]), 3)
        self.assertEqual(len(recipe.add_fragment_3p), 1)
        self.assertEqual(len(recipe.add_fragment_3p[0]), 4)

    def test_no_overhan(self):
        a = read_alignment(NO_OVERHANG_ALIGN)
        recipe = RecipeMaker(a).recipe
        self.assertEqual(len(recipe.add_fragment_5p), 0)
        self.assertEqual(len(recipe.add_fragment_3p), 0)


    def test_close_gaps(self):
        """Gaps close to each other should also work."""
        a = read_alignment(CLOSEGAPS)
        recipe = RecipeMaker(a).recipe
        self.assertEqual(len(recipe.difficult), 0)
        self.assertEqual(len(recipe.add_fragment), 4)


    def test_short_with_insert(self):
        a = read_alignment(WITH_INSERT)
        recipe = RecipeMaker(a).recipe
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
            recipe = RecipeMaker(align).recipe
            for ap, m in zip(align, matches):
                self.assertTrue(self.find_in_box(recipe, ap, m))



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


if __name__ == '__main__':
    main()

