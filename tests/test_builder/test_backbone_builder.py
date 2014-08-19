#!/usr/bin/env python
#
# test_backbone_builder.py
#
# unit tests for reconstructing full backbones
# e.g. after loop insertsions.
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
from moderna.ModernaStructure import ModernaStructure
from moderna.builder.BackboneBuilder import BackboneBuilder
from moderna.analyze.GeometryParameters import GeometryStandards
from test_data import *

DISCONNECTED = TEST_DATA_PATH + 'rna_structures/disconnected.pdb'
BB_MESSED_UP = TEST_DATA_PATH + 'rna_structures/bb_messed_up.pdb'

class BackboneBuilderTests(TestCase):
    def setUp(self):
        """Loads the A residue to start with."""
        self.struc = ModernaStructure('file',MINI_TEMPLATE)
        self.resi1 = self.struc['3']
        self.resi2 = self.struc['4']
        # delete the atoms that should be re-built.
        self.resi2.detach_child("O5'")
        self.resi2.detach_child("P")
        
    def tearDown(self):
        self.struc = None
        self.resi1 = None
        self.resi2 = None
        
    def test_validation(self):
        """Makes sure test data is set up properly."""
        self.assertFalse(self.resi2.is_backbone_intact())
        self.assertFalse(self.struc.are_residues_connected(self.resi1,self.resi2))
                         
    def test_build_backbone(self):
        """Checks whether the P+O5' atoms are constructed."""
        bb = BackboneBuilder(self.resi1,self.resi2, self.struc)
        self.assertTrue(self.resi2.child_dict.get('P'))
        self.assertTrue(self.resi2.child_dict.get("O5'"))

    def test_build_backbone_atomnames(self):
        """Checks whether the P+O5' atoms are constructed."""
        bb = BackboneBuilder(self.resi1,self.resi2, self.struc)
        self.assertTrue(self.resi2['P'])
        self.assertTrue(self.resi2["O5'"])
        self.assertEqual(self.resi2["P"].fullname, ' P')
        self.assertEqual(self.resi2["O5'"].fullname, " O5'")
        self.assertEqual(self.resi2["P"].element, 'P')
        self.assertEqual(self.resi2["O5'"].element, "O")
        

    def test_build_op1op2(self):
        """Checks whether the OP1 and OP2 atoms are constructed."""
        bb = BackboneBuilder(self.resi1,self.resi2, self.struc)
        self.assertTrue(self.resi2.child_dict.get('OP1'))
        self.assertTrue(self.resi2.child_dict.get("OP2"))

    def test_build_c5o5_distance(self):
        """Checks whether the C5'-O5' distance is OK."""
        bb = BackboneBuilder(self.resi1,self.resi2, self.struc)
        expected,stddev = GeometryStandards.bonds["X:O5',X:C5'"][0]
        self.assertTrue((self.resi2["C5'"]-self.resi2["O5'"]) - expected < stddev)

    def test_build_o5p_distance(self):
        """Checks whether the O5'-P distance is OK."""
        bb = BackboneBuilder(self.resi1,self.resi2, self.struc)
        expected,stddev = GeometryStandards.bonds["X:P,X:O5'"][0]
        self.assertTrue((self.resi2["O5'"]-self.resi2["P"]) - expected < stddev)

    def test_build_po3_distance(self):
        """Checks whether the P-O3' distance is OK."""
        bb = BackboneBuilder(self.resi1,self.resi2, self.struc)
        expected,stddev = GeometryStandards.bonds["X:O3',X+1:P"][0]
        self.assertTrue((self.resi2["P"]-self.resi1["O3'"]) - expected < stddev)
        
    def test_disconnected(self):
        """Structure example should not be OK"""
        # trna Leu by Marcin Skorupskiego
        struc = ModernaStructure('file',DISCONNECTED)
        self.assertFalse(struc.is_chain_continuous())
        bb = BackboneBuilder(struc['16'], struc['17'], struc)
        self.assertTrue(struc.is_chain_continuous())
        
    def _test_repair_congested(self):
        """Should avoid clashes between atoms."""
        #TODO: reactivate test
        struc = ModernaStructure('file',BB_MESSED_UP)
        self.assertFalse(struc['33'].is_phosphate_intact())
        bb = BackboneBuilder(struc['32'], struc['33'], struc)
        struc.write_pdb_file('out.pdb')
        self.assertTrue(bb.is_intact())
        #TODO: criteria are now more strict, therefore test fails.
        self.assertTrue(struc['32'].is_backbone_intact())
        self.assertTrue(struc['33'].is_backbone_intact())
        self.assertTrue(struc['33'].is_phosphate_intact())
        self.assertFalse(struc['33'].is_backbone_congested())
        self.assertFalse(struc['32'].is_backbone_congested())
        

    def test_build_conserve_c2c3(self):
        """downstream C2' and C3' positions should not change"""
        struc = ModernaStructure('file',FIXABLE_BACKBONE)
        resi1 = struc['4']
        resi2 = struc['5']
        c3_before = resi2["C3'"].get_vector()
        c2_before = resi2["C2'"].get_vector()
        bb = BackboneBuilder(resi1, resi2, struc)
        resi1 = struc['4']
        resi2 = struc['5']
        c3_after = resi2["C3'"].get_vector()
        c2_after = resi2["C2'"].get_vector()
        c2dist = c2_before-c2_after
        c3dist = c3_before-c3_after
        self.assertAlmostEqual(c2dist.norm(), 0.0, 3)
        self.assertAlmostEqual(c3dist.norm(), 0.0, 3)

    def test_build(self):
        struc = ModernaStructure('file',FIXABLE_BACKBONE)
        resi1 = struc['4']
        resi2 = struc['5']
        bb = BackboneBuilder(resi1, resi2, struc)
        #struc.write_pdb_file('out.pdb')
        self.assertTrue(bb.is_intact())
        
    def test_fccd(self):
        struc = ModernaStructure('file',FCCD_EXAMPLE)
        resi1 = struc['5']
        resi2 = struc['6']
        struc.remove_residue('5A')
        bb = BackboneBuilder(resi1, resi2, struc)
        self.assertTrue(bb.is_intact())
        
    def test_helix(self):
        struc = ModernaStructure('file', BROKEN_HELIX_EXAMPLE)
        resi1 = struc['0Z']
        resi2 = struc['1']
        bb = BackboneBuilder(resi1, resi2, struc)
        struc.write_pdb_file('out.pdb')
        self.assertTrue(bb.is_intact())
        
    def test_bad_insertion(self):
        """Should fix at least some breaks in a disastrous case"""
        struc = ModernaStructure('file', BAD_LOOP_INSERTION)
        breaks_before = str(struc.get_sequence()).count('_')
        self.assertEqual(breaks_before, 3)
        prev = None
        for resi in struc:
            if prev:
                bb = BackboneBuilder(prev, resi, struc)
            prev = resi
        breaks_after = str(struc.get_sequence()).count('_')
        self.assertTrue(breaks_after < breaks_before)
        
    

if __name__ == '__main__':
    main()
  
