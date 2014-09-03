#!/usr/bin/env python
#
# test_secstruc_calculator.py
#
# unit tests for calculating 2D structure from 3D
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

from unittest import main,TestCase
from moderna.ModernaStructure import ModernaStructure
from moderna.modifications import exchange_base
from test_data import *

class SecStrucCalculatorTests(TestCase):
    """Tests functions to calculate dot-bracket code."""
    def test_hairpin(self):
        """Calculates all pairs in a hairpin"""
        ms = ModernaStructure('file', RNA_HAIRPIN, 'D')
        ss = ms.get_secstruc()
        self.assertEqual(ss, '..((.......))')
        
    def test_1ehz(self):
        """Calculates all pairs in a tRNA"""
        # take note that 4-69 is a Wobble pair, and 19-56 would make a pseudoknot.
        expected = '(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))....'
        ms = ModernaStructure('file', RNA_1EHZ)
        ss = ms.get_secstruc()
        self.assertEqual(ss, expected)
        
    def test_pseudoknot(self):
        pass
        
    def test_replace(self):
        """Calculates secstruc in an edited structure."""
        ms = ModernaStructure('file', RNA_HAIRPIN, 'D')
        ss = ms.get_secstruc()
        self.assertEqual(ss, '..((.......))')
        exchange_base(ms['39'], 'G')
        ss = ms.get_secstruc()
        self.assertEqual(ss, '..((.......))')
        
    def test_nasty(self):
        """A secstruc is calculated from a trashy PDB."""
        ms = ModernaStructure('file', NASTY_PDB, 'A')
        ss = ms.get_secstruc()
        expected = '..(.((...((((........)))).((((.........)))).....(((((.......))))).)).)...............................................................................................................................................................................'
        self.assertEqual(ss, expected)
        

if __name__ == '__main__':
    main()
