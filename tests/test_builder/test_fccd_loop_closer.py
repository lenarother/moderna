#!/usr/bin/env python
"""
Tests for the FCCDLoopCloser module.
"""

from unittest import TestCase, main
from moderna.builder.FCCDLoopCloser import FCCDLoopCloser
from moderna.tests.test_data import FCCD_EXAMPLE
from Bio.PDB import PDBParser
from Bio.PDB.Vector import Vector
import numpy


class FCCDLoopCloserTests(TestCase):
    """
    Tests the FCCD loop closing algorithm to model gaps between two riboses.
    """
    def setUp(self):
        """Initialize test data."""
        parser = PDBParser()
        struc = parser.get_structure('riboses',FCCD_EXAMPLE)
        resi5  = struc[0]['A'][(' ',5,' ')]
        resi5a = struc[0]['A'][(' ',5,'A')]
        resi6  = struc[0]['A'][(' ',6,' ')]

        self.fixed = [resi6["C4'"], resi6["C3'"], resi6["C2'"]]
        self.moving = [resi5["O4'"], resi5["C4'"], \
                       resi5["C3'"], resi5["O3'"], \
                       resi5a["P"], resi5a["O5'"], resi5a["C5'"], \
                       resi5a["C4'"], resi5a["C3'"], resi5a["C2'"]
                       ]
        self.fccd = FCCDLoopCloser(self.moving, self.fixed)
        
    def tearDown(self):
        """Cleans up test data."""
        self.rna = None
        self.fixed = None
        self.moving = None
        self.fccd = None
        
    def test_small_rmsd(self):
        """After running FCCD, the example should have a small RMSD."""
        msg, rmsd, n_it = self.fccd.run_fccd(threshold=0.1)
        self.fccd.copy_vectors_to_atoms()
        self.assertTrue(msg.startswith('RMSD threshold reached'))
        self.assertTrue(rmsd < 0.1)
        self.assertTrue(n_it > 0)
        
    def test_close_loop(self):
        """After running RMSD, the gap in the riboses should be closed."""
        #TODO: uncomment after contributing RNA structure code.
        #self.assertFalse(self.rna.are_residues_connected(\
        #        self.rna['5'], self.rna['6']))
        self.fccd.run_fccd(threshold=0.1)
        self.fccd.copy_vectors_to_atoms()
        #self.assertTrue(self.rna.are_residues_connected(\
        #        self.rna['5'], self.rna['6']))

    def test_get_moving_coords(self):
        """Method for moving coordinates should pass."""
        center = Vector([100.0, 0.0, 0.0])
        self.fccd.get_moving_coords(center)
        self.assertTrue(True)
        self.assertTrue(True)

    def test_get_fixed_coords(self):
        """Method for fixed coordinates should pass."""
        center = Vector([100.0, 0.0, 0.0])
        self.fccd.get_fixed_coords(center)
        self.assertTrue(True)
        
    def test_copy_vectors(self):
        """Should modify the atom coordinates."""
        self.fccd.moving[0] = Vector([1.23, 4.56, 7.89])
        self.fccd.copy_vectors_to_atoms()
        self.assertEqual(self.moving[0].coord[0], 1.23)
        self.assertEqual(self.moving[0].coord[1], 4.56)
        self.assertEqual(self.moving[0].coord[2], 7.89)

if __name__ == '__main__':
    main()
