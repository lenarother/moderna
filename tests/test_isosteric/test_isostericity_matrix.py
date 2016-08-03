#!/usr/bin/env python
"""
unit tests for different functions calculating isosteric base pairs
"""

from unittest import main, TestCase
from moderna.isosteric.IsostericityMatrices import IsostericityMatrices


class IsostericityMatrixTests(TestCase):
    
    def test_check_isostericity(self):
        """aaa"""
        im = IsostericityMatrices()
        result = im.check_isostericity('AC','GU','cWW')
        result2 = im.check_isostericity('AA','GG','cWW')
        result3 = im.check_isostericity('AC','UU','cWW')
        self.assertTrue(result)
        self.assertFalse(result2)
        self.assertFalse(result3)

    def test_check_isostericity_cutoff(self):
        """aaa"""
        im = IsostericityMatrices()
        result = im.check_isostericity('UU','AA','tWH',  2.0)
        result2 = im.check_isostericity('UU','AA','tWH',  3.0)
        self.assertFalse(result)
        self.assertTrue(result2)
        
    def test_check_isostericity_none(self):
        """aaa"""
        im = IsostericityMatrices()
        result = im.check_isostericity('AA','GA','tWS')
        self.assertFalse(result)

    def test_show_isosteric_bp(self):
        """aaa"""
        im = IsostericityMatrices()
        result = im.show_isosteric_bp('AC','cWW')
        self.assertEqual(result, ('AC','GU'))
   
    def test_show_isosteric_bp_cutoff(self):
       """aaa"""
       im = IsostericityMatrices()
       result = im.show_isosteric_bp('AC', 'cWW', 3.0)
       result = list(result)
       result.sort()
       self.assertEqual(result, ['AC', 'AU', 'CG', 'GC', 'GU', 'UA', 'UU'])
   
if __name__ == '__main__':
    main()
    
