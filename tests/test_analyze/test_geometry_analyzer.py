#!/usr/bin/env python
"""
Unit Tests for GeometryAnalyzer
"""


from unittest import main, TestCase
from moderna.analyze.GeometryAnalyzer import GeometryAnalyzer
from moderna.ModernaStructure import ModernaStructure
from moderna.commands import create_model
from test_data import MINI_TEMPLATE, BAD_TEMPLATE


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
        self.assertTrue('Unusual bond lengths:' not in r)
        self.assertTrue('Unusual bond angles:' not in r)
        self.assertTrue('Unusual dihedral angles:' in r)
        self.bad_ga.analyze()
        r = str(self.bad_ga)
        self.assertTrue('Unusual bond lengths:' in r)
        self.assertTrue('Unusual bond angles:' in r)
        # check explicit values
        self.assertTrue('2.58' in r)
        self.assertTrue('151.85' in r)
        self.assertTrue('248.03' in r)

    def test_empty(self):
        """An empty structure creates an empty report."""
        m = create_model()
        ga = GeometryAnalyzer(m)
        self.assertEqual("", str(ga).strip())

if __name__ == '__main__':
    main()
    
