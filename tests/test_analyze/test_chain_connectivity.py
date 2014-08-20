
from unittest import main, TestCase
from Bio.PDB.PDBParser import PDBParser
from moderna.analyze.ChainConnectivity import are_residues_connected, \
    is_backbone_complete, is_backbone_intact, \
    is_phosphate_intact, is_backbone_congested
from moderna.ModernaStructure import ModernaStructure
from moderna.RNAResidue import RNAResidue
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


class ResidueIntegrityTests(TestCase):

    def setUp(self):
        """Loads the A residue to start with."""
        self.a=PDBParser().get_structure('test_struc',A_RESIDUE)[0].child_list[0].child_list[0]
        self.chain=PDBParser().get_structure('test_struc',MINI_TEMPLATE)[0].child_list[0]

    def test_is_backbone_complete(self):
        """Complete backbone atoms are recognized."""
        for resi in self.chain:
            resi = RNAResidue(resi)
            self.assertTrue(is_backbone_complete(resi))

    def test_backbone_incomplete(self):
        """Negative example with missing backbone atoms."""
        chain=PDBParser().get_structure('test_struc',INCOMPLETE_BACKBONE)[0].child_list[0]
        for resi in chain:
            resi = RNAResidue(resi)
            self.assertFalse(is_backbone_complete(resi))

    def test_is_backbone_intact(self):
        """Check all kinds of backbone discontinuities in one residue."""
        chain=PDBParser().get_structure('test_struc',BROKEN_BACKBONE)[0].child_list[0]
        residues = [r for r in chain]
        for resi in residues[:5]:
            mr = RNAResidue(resi)
            self.assertFalse(is_backbone_intact(mr))
        mr = RNAResidue(chain[6])
        self.assertTrue(is_backbone_intact(mr))

    def test_is_backbone_intact_5p3p(self):
        """Check all kinds of backbone discontinuities in one residue."""
        chain=PDBParser().get_structure('test_struc',BROKEN_BACKBONE)[0].child_list[0]
        residues = [r for r in chain]
        result_5p = []
        result_3p = []
        for resi in residues[:6]:
            mr = RNAResidue(resi)
            result_5p.append(is_backbone_intact(mr, mode="5'"))
            result_3p.append(is_backbone_intact(mr, mode="3'"))
        self.assertEqual(result_5p, [False, False, False, False, True, True])
        self.assertEqual(result_3p, [True, True, True, False, False, True])

    def test_is_phosphate_intact(self):
        """Check whether OP1 and OP2 are in place"""
        chain=PDBParser().get_structure('test_struc',BB_MESSED_UP)[0].child_list[0]
        resi1 = RNAResidue(chain[('H_c  ', 32, ' ')])
        resi2 = RNAResidue(chain[(' ', 33, ' ')])
        self.assertTrue(is_phosphate_intact(resi1))
        self.assertFalse(is_phosphate_intact(resi2))

    def test_is_backbone_congested(self):
        """Check whether backbone atoms clash into rest of the structure."""
        resi = RNAResidue(self.chain.child_list[2])
        self.assertFalse(is_backbone_congested(resi))
        # now check a structure where the backbone clashes into O2'
        chain=PDBParser().get_structure('test_struc', BB_MESSED_UP)[0].child_list[0]
        resi = RNAResidue(chain[('H_c  ', 32, ' ')])
        self.assertTrue(is_backbone_congested(resi))
        resi = RNAResidue(chain[(' ', 33, ' ')])
        self.assertTrue(is_backbone_congested(resi))


if __name__ == '__main__':
    main()

