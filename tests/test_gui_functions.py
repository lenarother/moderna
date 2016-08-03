#!/usr/bin/env python
"""
Unit Tests for helper functions for the ModernaGUI
"""

from unittest import main, TestCase
from test_data import *
from moderna import *
from moderna.analyze.topology_matcher import topology_graph
from numpy import array, ndarray


class BondCoordTests(TestCase):
    """Tests for the function delivering coord tuples for bonds."""

    def setUp(self):
        self.t = Template(MINI_TEMPLATE)
        self.a = Template(A_RESIDUE)
        
    def test_get_bond_coordinates(self):
        """Should create a list of BondTopology objects."""
        result = topology_graph.get_bond_coordinates(self.t)
        self.assertTrue(isinstance(result, list))
        self.assertTrue(len(result) > 0 )
        
    def test_get_bond_coord_singleresi(self):
        result = topology_graph.get_bond_coordinates(self.a)
        self.assertEqual(len(result), 24)

    def test_get_bond_coord_struc(self):
        result = topology_graph.get_bond_coordinates(self.t)
        self.assertEqual(len(result), 360)

    def test_get_bond_coord_with_coordinates(self):
        """Should create numpy.arrays in result"""
        result = topology_graph.get_bond_coordinates(self.a)
        first = result[0]
        self.assertTrue(isinstance(first, topology_graph.BondTopology))
        self.assertTrue(isinstance(first[0], ndarray))
        self.assertTrue(isinstance(first[1], ndarray))
        

if __name__ == '__main__':
    main()
    
