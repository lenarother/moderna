#!/usr/bin/env python
#
# test_geometry_analyzer.py
#
# unit tests GeometryAnalyzer
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__contributors__ = "Magdalena Rother, Kristian Rother"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"


from unittest import main, TestCase
from moderna.analyze.GeometryAnalyzer import GeometryAnalyzer
from moderna.ModernaStructure import ModernaStructure
from moderna.commands import create_model
from test_data import *


class GeometryAnalyzerTests(TestCase):

    def setUp(self):
        t = ModernaStructure('file', MINI_TEMPLATE)
        self.ga = GeometryAnalyzer(t)
        t = ModernaStructure('file', BAD_TEMPLATE, chain_name='A', seq=t.get_sequence())
        self.bad_ga = GeometryAnalyzer(t)
        
    def test_check_bonds(self):
        """Should return a list of bad bond lengths."""
        r = self.ga.check_bonds()
        self.assertEqual(len(r), 0)
        r = self.bad_ga.check_bonds()
        self.assertEqual(len(r), 3)
        
    def test_check_angles(self):
        """Should return a list of bad angles."""
        r = self.ga.check_angles()
        self.assertEqual(len(r), 0)
        r = self.bad_ga.check_angles()
        self.assertEqual(len(r), 5)
        
    def test_check_dihedrals(self):
        """Should return a list of bad dihedrals."""
        r = self.ga.check_dihedrals()
        self.assertEqual(len(r), 1)
        r = self.bad_ga.check_dihedrals()
        self.assertEqual(len(r), 2)
    
    def test_report(self):
        """The report should only contain sections with a notification."""
        self.ga.analyze()
        r = str(self.ga)
        self.assertTrue(r.find('Unusual bond lengths:')==-1)
        self.assertTrue(r.find('Unusual bond angles:')==-1)
        self.assertTrue(r.find('Unusual dihedral angles:')>-1)
        self.bad_ga.analyze()
        r = str(self.bad_ga)
        self.assertTrue(r.find('Unusual bond lengths:')>-1)
        self.assertTrue(r.find('Unusual bond angles:')>-1)
        # check explicit values
        self.assertTrue(r.find('2.58')>-1)
        self.assertTrue(r.find('151.85')>-1)
        self.assertTrue(r.find('248.03')>-1)
        
    def test_empty(self):
        """An empty structure creates an empty report."""
        m = create_model()
        ga = GeometryAnalyzer(m)
        self.assertEqual("", str(ga).strip())

if __name__ == '__main__':
    main()
    
