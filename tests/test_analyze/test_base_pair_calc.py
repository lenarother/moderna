#!/usr/bin/env python
"""
unit tests for different functions calculating base pairs.
"""

from unittest import main, TestCase
from moderna.ModernaStructure import  ModernaStructure
from test_data import *

EXAMPLE_CHAR =  TEST_DATA_PATH+'rna_structures/rna_with_char_resid_1n77_C.pdb'


class BasePairCalcTests(TestCase):
    
    def test_wc_pairs(self):
        """Should recognize W-C-pairs"""
        m = ModernaStructure('file', RNA_HAIRPIN, 'D')
        self.assertEqual(m['31'].get_bp(m['39']).type, '+/+')
        self.assertEqual(m['39'].get_bp(m['31']).type, '+/+')
        self.assertEqual(m['30'].get_bp(m['40']).type, '+/+')
        self.assertEqual(m['40'].get_bp(m['30']).type, '+/+')
        
    def test_residue_no_pair(self):
        """Two unpaired residues result in None"""
        m = ModernaStructure('file', RNA_HAIRPIN, 'D')
        result = m['30'].get_bp(m['38'])
        self.assertEqual(result, None)
        
    def test_get_base_pairs(self):
        """Calculates one canonical and two noncanonicals correctly."""
        m = ModernaStructure('file', RNA_HAIRPIN, 'D')
        result = m.get_base_pairs()
        result = map(str, result['31'])
        self.assertTrue('31 +/+ 39' in result)  
        self.assertTrue('31 HS 32' in result)  
        self.assertTrue('31 SS 30' in result)  
        
    def test_structure_unpaired(self):
        """Single-stranded structure has no pairs."""
        m = ModernaStructure('file', EXAMPLE_CHAR, 'C')
        result = m.get_base_pairs()
        self.assertEqual(result['37'], [])                

    def test_get_base_pairs_generator(self):
        """aaa"""
        m = ModernaStructure('file', RNA_HAIRPIN, 'D')
        result = [str(x) for x in m.get_base_pairs('generator')]
        self.assertTrue('30 +/+ 40' in result)
        self.assertTrue(len(result) >= 4)
        # TODO: check with RNAView output.
        
    def test_wc_pairs_hairpin(self):
        """Exactly two W-C pairs in a hairpin are recognized."""
        m = ModernaStructure('file', RNA_HAIRPIN, 'D')
        result = m.get_base_pairs()
        # calc number of +/+ pairs
        count = 0
        for val in result.values():
            for bpair in val:
                if bpair.type == '+/+':
                    count += 1
        self.assertEqual(count, 4)
    
    def test_wc_pairs_1ehz(self):
        """All W-C pairs in a tRNA are recognized."""
        m = ModernaStructure('file', RNA_1EHZ, 'A')
        result = m.get_base_pairs()
        # calc number of +/+ pairs
        gc_count, au_count = 0, 0
        for key in result:
            val = result[key]
            for bpair in val:
                if bpair.type == '+/+':
                    gc_count += 1
                elif bpair.type == '-/-':
                    au_count += 1
        #print gc_count, au_count
        self.assertEqual(au_count, 14)
        self.assertEqual(gc_count, 26)
        """# canonical from 1ehz according to rnaview
     72, A:     1 G-C    72 A: +/+ cis         XIX
     2_71, A:     2 C-G    71 A: +/+ cis         XIX
     3_70, A:     3 G-C    70 A: +/+ cis         XIX
     4_69, A:     4 G-U    69 A: W/W cis         XXVIII
     5_68, A:     5 A-U    68 A: -/- cis         XX
     6_67, A:     6 U-A    67 A: -/- cis         XX
     7_66, A:     7 U-A    66 A: -/- cis         XX
    10_25, A:    10 g-C    25 A: +/+ cis         XIX
    11_24, A:    11 C-G    24 A: +/+ cis         XIX
    12_23, A:    12 U-A    23 A: -/- cis         XX
    13_22, A:    13 C-G    22 A: +/+ cis         XIX
    19_56, A:    19 G-C    56 A: +/+ cis    syn    XIX
    27_43, A:    27 C-G    43 A: +/+ cis         XIX
    28_42, A:    28 C-G    42 A: +/+ cis         XIX
    29_41, A:    29 A-U    41 A: -/- cis         XX
    30_40, A:    30 G-c    40 A: +/+ cis         XIX
    49_65, A:    49 c-G    65 A: +/+ cis         XIX
    50_64, A:    50 U-A    64 A: -/- cis         XX
    51_63, A:    51 G-C    63 A: +/+ cis         XIX
    52_62, A:    52 U-A    62 A: -/- cis         XX
    53_61, A:    53 G-C    61 A: +/+ cis         XIX
    """

if __name__ == '__main__':
    main()
    
