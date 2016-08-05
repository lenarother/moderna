#!/usr/bin/env python
"""
Unit Tests for the StackingCalculator class
"""

from unittest import TestCase, main
from moderna.analyze.StackingCalculator import StackingCalculator
from moderna.ModernaStructure import ModernaStructure
from moderna.tests.test_data import *

EXPECTED_STACKING = [
    ('1','2','>>','<<'),
    ('3','4','>>','<<'),
    ('4','5','>>','<<'),
    ('5','6','>>','<<'),
    ('6','7','>>','<<'),
    ('8','13','><','><'),
    ('10','11','>>','<<'),
    ('11','12','>>','<<'),
    ('12','13','>>','<<'),
    ]
    
class StackingCalculatorTests(TestCase):
    
    def setUp(self):
        self.s = ModernaStructure('file',MINI_TEMPLATE)
        self.sc = StackingCalculator()

    def test_get_stacking(self):
        """Returned stacking values should correspond to those given."""
        for stacking in self.sc.get_stacking(self.s):
            found = False
            for ref1, ref2, stack, revstack in EXPECTED_STACKING:
                if (stacking.resi1.identifier == ref1 and stacking.resi2.identifier == ref2):
                    self.assertEqual(stacking.type,stack)
                    found = True
                elif (stacking.resi2.identifier == ref1 and stacking.resi1.identifier == ref2):
                    self.assertEqual(stacking.type,revstack)
                    found = True
                    
            self.assertTrue(found)
            
    def test_get_stacking_ehz(self):
        m = ModernaStructure('file', RNA_1EHZ, 'A')
        result = self.sc.get_stacking(m)
        self.assertEqual(len(result), 57)


if __name__ == '__main__':
    main()
