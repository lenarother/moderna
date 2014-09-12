#!/usr/bin/env python
#
# test_chi_rotation.py
#
# unit tests for rotating base around the glycosidic bond
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

from unittest import TestCase, main
from moderna.ModernaStructure import ModernaStructure
from moderna.builder.ChiRotator import rotate_chi
from test_data import *
from Bio.PDB.Vector import calc_dihedral
import math

class ChiRotationTests(TestCase):

    def setUp(self):
        self.s = ModernaStructure('file',MINI_TEMPLATE)
        self.r = self.s['3']
        self.mod = self.s['10']
        self.ehz = ModernaStructure('file',RNA_1EHZ)

    def calc_chi(self,r):
        if r.long_abbrev in ['C','U']:
            n1,n2 = ('N1','C2')
        else:
            n1,n2 = ('N9','C8')
        chi = calc_dihedral(r["C2'"].get_vector(), r["C1'"].get_vector(),\
            r[n1].get_vector(),r[n2].get_vector())
        return chi * 180/math.pi

    def test_rotate_dihedral(self):
        """Dihedral angle of glycosidic bonds should change."""
        chi_before = self.calc_chi(self.r)
        rotate_chi(self.r, 10.0)
        chi_after = self.calc_chi(self.r)
        chi_diff = chi_after-chi_before
        self.assertAlmostEqual(chi_diff,10.0,2)

    def test_rotate_base(self):
        """Atoms of the base should change."""
        changing = ["N1","C2","N2","N3","C4","C5","C6","N7","C8"]
        for aname in changing:
            coord_before = self.r[aname].coord
            rotate_chi(self.r, 5.0)
            coord_after = self.r[aname].coord
            self.assertNotEqual(list(coord_after), list(coord_before))

    def test_not_rotate_sugar(self):
        """Atoms of the ribose plus glycosidic bond should not change."""
        not_changing = ["C1'","N9","C2'","O2'","C3'","C4'","C5'","O3'","O4'","O5'"]
        for aname in not_changing:
            coord_before = self.r[aname].coord
            rotate_chi(self.r, 5.0)
            coord_after = self.r[aname].coord
            self.assertEqual(list(coord_after), list(coord_before))

    def test_rotate_modifications(self):
        """Modified atoms on the base should be rotated as well."""
        coord_before = self.mod['CM2'].coord
        rotate_chi(self.mod, 50.0)
        coord_after = self.mod['CM2'].coord
        self.assertNotEqual(list(coord_after), list(coord_before))

    def test_not_rotate_methyl(self):
        """Modifications of the O2' should not be rotated."""
        mod = self.ehz['32']
        coord_before = mod['CM2'].coord
        rotate_chi(self.mod, 50.0)
        coord_after = mod['CM2'].coord
        self.assertEqual(list(coord_after), list(coord_before))
        

if __name__ == '__main__':
    main()
