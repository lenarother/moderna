#!/usr/bin/env python
#
# test_pucker_calculator.py
#
# unit tests for pucker calculation
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

from unittest import TestCase,main
from moderna.ModernaStructure import ModernaStructure
from moderna.analyze.PuckerCalculator import PuckerCalculator
from test_data import RNA_1EHZ

PUCKER_EXAMPLES = [
    ('11',"C3'-endo"),
    ('12',"C3'-endo"),
    ('13',"C3'-endo"),
    ('14',"C2'-exo"),
    ('16',"C3'-endo"),
    ]

class PuckerCalculatorTests(TestCase):
    """Tests for calculating ribose puckers."""
    def setUp(self):
        """Initializes class instances used for testing."""
        self.struc = ModernaStructure('file', RNA_1EHZ)

    def test_get_pucker_for_residue(self):
        """Should return a string with sugar pucker."""
        for resi,pucker in PUCKER_EXAMPLES:
            pc = PuckerCalculator()
            puk = pc.get_pucker(self.struc[resi])
            self.assertEqual(puk, pucker)

    def test_pucker_property(self):
        """Pucker is a property of ModernaResidue."""
        for resi,pucker in PUCKER_EXAMPLES:
            self.assertEqual(self.struc[resi].pucker, pucker)

if __name__ == '__main__':
    main()
    
