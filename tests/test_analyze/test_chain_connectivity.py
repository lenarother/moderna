
from unittest import main, TestCase
from moderna.analyze.ChainConnectivity import are_residues_connected
from moderna.ModernaStructure import ModernaStructure
from moderna.sequence.ModernaSequence import Sequence
from test_data import *


class ChainConnectivityTests(TestCase):

    def setUp(self):
        """Loads a structure to start with."""
        self.t = ModernaStructure('file',MINI_TEMPLATE)

    def test_is_connected_to_true(self):
        """The 3'O (n) and P (n+1) must be close together."""
        connected = are_residues_connected(self.t['4'], self.t['5'])
        self.assertTrue(connected)

    def test_is_connected_to_false(self):
        """The 3'O (n) and P (n+1) must be close together."""
        connected = are_residues_connected(self.t['4'], self.t['6'])
        self.assertFalse(connected)

    def test_is_connected_reverse(self):
        """Reverse order of residues changes the result."""
        connected = are_residues_connected(self.t['5'], self.t['4'])
        self.assertFalse(connected)




if __name__ == '__main__':
    main()

