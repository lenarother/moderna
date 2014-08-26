#!/usr/bin/env python
#
# test_backbone_builder.py
#
# unit tests for reconstructing P and O5' atoms
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
from moderna.builder.PhosphateBuilder import PhosphateBuilder
from moderna.analyze.GeometryParameters import GeometryStandards
from moderna.analyze.ChainConnectivity import are_residues_connected, is_backbone_intact
from Bio.PDB.Vector import Vector, calc_angle, calc_dihedral
import math
from test_data import *


class PhosphateBuilderTests(TestCase):
    def setUp(self):
        """Loads the A residue to start with."""
        self.struc = ModernaStructure('file',MINI_TEMPLATE)
        self.resi1 = self.struc['3']
        self.resi2 = self.struc['4']
        # delete the atoms that should be re-built.
        self.resi2.detach_child("O5'")
        self.resi2.detach_child("P")
        
    def test_validation(self):
        """Makes sure test data is set up properly."""
        self.assertFalse(is_backbone_intact(self.resi2))
        self.assertFalse(are_residues_connected(self.resi1,self.resi2))
                         
    def test_build_backbone(self):
        """Checks whether the P+O5' atoms are constructed."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        self.assertTrue(self.resi2.child_dict.get('P'))
        self.assertTrue(self.resi2.child_dict.get("O5'"))

    def test_build_backbone_atomnames(self):
        """Checks whether the P+O5' atoms are constructed."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        self.assertTrue(self.resi2['P'])
        self.assertTrue(self.resi2["O5'"])
        self.assertEqual(self.resi2["P"].fullname, ' P')
        self.assertEqual(self.resi2["O5'"].fullname, " O5'")
        self.assertEqual(self.resi2["P"].element, 'P')
        self.assertEqual(self.resi2["O5'"].element, "O")
        

    def test_build_op1op2(self):
        """Checks whether the OP1 and OP2 atoms are constructed."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        self.assertTrue(self.resi2.child_dict.get('OP1'))
        self.assertTrue(self.resi2.child_dict.get("OP2"))

    def test_build_c5o5_distance(self):
        """Checks whether the C5'-O5' distance is OK."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        low,high = GeometryStandards.bonds["X:O5',X:C5'"][0]
        self.assertTrue(low <= (self.resi2["C5'"]-self.resi2["O5'"]) <= high)

    def test_build_o5p_distance(self):
        """Checks whether the O5'-P distance is OK."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        low, high = GeometryStandards.bonds["X:P,X:O5'"][0]
        self.assertTrue(low <= (self.resi2["O5'"]-self.resi2["P"]) <= high)

    def test_build_po3_distance(self):
        """Checks whether the P-O3' distance is OK."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        low, high = GeometryStandards.bonds["X:O3',X+1:P"][0]
        self.assertTrue(low <= (self.resi2["P"]-self.resi1["O3'"]) <= high)

    def test_c3o3p_angle(self):
        """Checks whether the C3'-O3'-P angle is OK."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        amin,amax = GeometryStandards.angles["X:C3',X:O3',X+1:P"][0]
        angle = calc_angle(self.resi1["C3'"].get_vector(),self.resi1["O3'"].get_vector(),self.resi2["P"].get_vector())
        self.assertTrue(amin < math.degrees(angle) < amax)

    def test_o3po5_angle(self):
        """Checks whether the O3'-P-O5' angle is OK."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        amin,amax = GeometryStandards.angles["X:O3',X+1:P,X+1:O5'"][0]
        angle = calc_angle(self.resi1["O3'"].get_vector(),self.resi2["P"].get_vector(),self.resi2["O5'"].get_vector())
        self.assertTrue(amin < math.degrees(angle) < amax)

    def test_po5c5_angle(self):
        """Checks whether the P-O5'-C5' angle is OK."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        amin,amax = GeometryStandards.angles["X:P,X:O5',X:C5'"][0]
        angle = calc_angle(self.resi2["P"].get_vector(),self.resi2["O5'"].get_vector(),self.resi2["C5'"].get_vector())
        self.assertTrue(amin < math.degrees(angle) < amax)

    def test_o5c5c4_angle(self):
        """Checks whether the O5'-C5'-C4' angle is OK."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        amin,amax = GeometryStandards.angles["X:O5',X:C5',X:C4'"][0]
        angle = calc_angle(self.resi2["O5'"].get_vector(),self.resi2["C5'"].get_vector(),self.resi2["C4'"].get_vector())
        self.assertTrue(amin < math.degrees(angle) < amax)

    def test_build_backbone_connection(self):
        """Checks whether two residues are properly connected."""
        bb = PhosphateBuilder(self.resi1,self.resi2)    
        bb.build()
        self.assertTrue(are_residues_connected(self.resi1,self.resi2))
        
    def test_op1op2_geometry(self):
        """Checks whether the OP1 and OP2 geometry is OK."""
        bb = PhosphateBuilder(self.resi1,self.resi2)
        bb.build()
        self.assertTrue(1.2 < self.resi2['P']-self.resi2['OP1'] < 1.8)
        self.assertTrue(1.2 < self.resi2['P']-self.resi2['OP2'] < 1.8)
        self.assertTrue(self.resi2['OP1']-self.resi2['OP2'] > 1.5)
        self.assertTrue(self.resi2['OP1']-self.resi2["O5'"] > 1.5)
        self.assertTrue(self.resi2['OP2']-self.resi2["O5'"] > 1.5)
        self.assertTrue(self.resi2['OP1']-self.resi1["O3'"] > 1.5)
        self.assertTrue(self.resi2['OP2']-self.resi1["O3'"] > 1.5)
'''
    def test_mm(self):
        m = ModernaStructure('file', 'm1.pdb')
        bb = BackboneBuilder(m['0C'], m['1'])
'''
if __name__ == '__main__':
    main()
  
