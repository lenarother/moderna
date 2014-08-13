#!/usr/bin/env python
#
# test_chain_discontinuity.py
#
# unit tests for recognizing connected and disconnected residues.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

from unittest import main, TestCase
from moderna.ModernaStructure import ModernaStructure
from moderna.ModernaResidue import ModernaResidue
from moderna.ModernaSequence import Sequence
from test_data import *


class ChainDiscontinuityTests(TestCase):
    def setUp(self):
        """Loads a structure to start with."""
        self.t = ModernaStructure('file',MINI_TEMPLATE)

    def test_is_connected_to_true(self):
        """The 3'O (n) and P (n+1) must be close together."""
        connected = self.t.are_residues_connected(self.t['4'],self.t['5'])
        self.assertTrue(connected)

    def test_is_connected_to_false(self):
        """The 3'O (n) and P (n+1) must be close together."""
        connected = self.t.are_residues_connected(self.t['4'],self.t['6'])
        self.assertFalse(connected)
        
    def test_is_connected_strict(self):
        """All backbones need to be connected between two bases."""
        seq = self.t.get_sequence()
        s = ModernaStructure('file',BROKEN_BACKBONE)
        self.assertEqual(s.get_sequence(), Sequence('G_C_GG_A_U_UUALCUCAG'))
        
    


if __name__ == '__main__':
    main()
  
