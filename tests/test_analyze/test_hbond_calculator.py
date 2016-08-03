#!/usr/bin/env python
"""
Unit Tests for checking interactions calculation
"""

from moderna.analyze.HBondCalculator import HBondCalculator
from unittest import main, TestCase
from moderna.ModernaStructure import ModernaStructure
from test_data import *
import os


class HBondCalculatorTests(TestCase):

    def setUp(self):
        """Initializes class instances used for testing."""
        self.rna = ModernaStructure('file',RNA_TWO_PIECES)
        self.example1 = ModernaStructure('file',AC_BASEPAIR, '0')
        self.example1_dist = ModernaStructure('file',DISTANT_AC_BASEPAIR, '0')
        self.ehz = ModernaStructure('file', RNA_1EHZ, 'A')
        self.hb = HBondCalculator()

    def test_calc_hbond_list(self):
        """Should calculate list of H-bonds from two residues."""
        res1,res2 = self.rna['1'],self.rna['72']
        hbonds = self.hb.calc_hbond_list(res1,res2)
        self.assertEqual(len(hbonds),3)

    def test_calc_hbond(self):
        """Should calculate H-bond from two atoms."""
        acceptor = self.rna['1'].child_dict["O6"]
        donor = self.rna['72'].child_dict["N4"]
        hbond = self.hb.calc_hbond(donor,acceptor)
        self.assertNotEqual(hbond,None)
        
    def test_calc_hbond_modified(self):
        """Calculates Hbonds of modified residues."""
        hbonds = self.hb.calc_hbond_list(self.ehz['25'], self.ehz['10'])
        self.assertNotEqual(hbonds,[])
        

    def test_calc_hbond_fail(self):
        """Should not calculate H-bond for atoms with strange angles."""
        atom1 = self.rna['1'].child_dict["N1"]
        atom2 = self.rna['72'].child_dict["N4"]
        hbond = self.hb.calc_hbond(atom2,atom1)
        self.assertEqual(hbond,None)

    def test_calc_hbond_dist(self):
        """Should not calculate H-bond for atoms with long distance."""
        atom1 = self.rna['1'].child_dict["O6"]
        atom2 = self.rna['70'].child_dict["N4"]
        hbond = self.hb.calc_hbond(atom2,atom1)
        self.assertEqual(hbond,None)
        
    def test_get_hbond_donors(self):
        """Should get all donor atoms."""
        expected = ["N1","C8","N2","O2'",]
        don = list(self.rna['1'].get_hbond_donors())
        self.assertEqual(len(don),len(expected))
        for d in don:
            self.assertTrue(d.name in expected)
            
    def test_get_hbond_acceptors(self):
        """Should get all acceptor atoms."""
        expected = ["O2'","O3'","O4'","O5'","OP1","OP2","N9","N7","O6","N3","N1"]
        acc = list(self.rna['1'].get_hbond_acceptors())
        self.assertEqual(len(acc),len(expected))
        for a in acc:
            self.assertTrue(a.name in expected)
            
    def test_calc_hbond_angle(self):
        """Hbond should depend on the angle"""
        don, acc = self.example1['498']['N6'], self.example1['494']['O2']
        self.assertTrue(self.hb.calc_hbond(don, acc))
        don, acc = self.example1_dist['498']['N6'], self.example1_dist['494']['O2']
        self.assertFalse(self.hb.calc_hbond(don, acc))
        
    def test_get_neighbors_purine(self):
        """Should return a list of close atoms"""
        #TODO: transfer to ModernaResidue
        a1 = self.example1['498']['N6']
        nb = self.example1['498'].get_neighbors(a1)
        self.assertEqual(len(nb), 1)
        self.assertTrue(self.example1['498']['C6'] in nb)
        # second example
        a1 = self.example1['498']['N9']
        nb = self.example1['498'].get_neighbors(a1)
        self.assertEqual(len(nb), 3)
        self.assertTrue(self.example1['498']["C1'"] in nb)
        self.assertTrue(self.example1['498']["C4"] in nb)
        self.assertTrue(self.example1['498']["C8"] in nb)

    def test_get_neighbors_pyrimidine(self):
        """Should return a list of close atoms"""
        #TODO: transfer to ModernaResidue
        a1 = self.rna['72']['C6']
        nb = self.rna['72'].get_neighbors(a1)
        self.assertEqual(len(nb), 2)
        self.assertTrue(self.rna['72']['N1'] in nb)
        self.assertTrue(self.rna['72']['C5'] in nb)
        # one more
        a1 = self.rna['72']['C4']
        nb = self.rna['72'].get_neighbors(a1)
        self.assertEqual(len(nb), 3)
        self.assertTrue(self.rna['72']['N3'] in nb)
        self.assertTrue(self.rna['72']['N4'] in nb)
        self.assertTrue(self.rna['72']['C5'] in nb)
        

    def test_get_hydrogens(self):
        """Should construct the right number of hydrogen positions"""
        #TODO: transfer to ModernaResidue        
        o2 = self.example1['498']["O2'"]
        hydro = [h for h in self.example1['498'].get_donor_hydrogens(o2)]
        self.assertEqual(len(hydro), 12)
        c8 = self.example1['498']["C8"]
        hydro = [h for h in self.example1['498'].get_donor_hydrogens(c8)]
        self.assertEqual(len(hydro), 1)


class HBondWesthofTests(TestCase):
 
    def setUp(self):
        """Initializes class instances used for testing."""
        self.hb = HBondCalculator()
        
    def check_example(self, fn):
            ll = open(fn).readlines()
            ll = [l for l in ll if l.startswith('ATOM')]
            chain1 = ll[3][21]
            chain2 = ll[-3][21]
            bp1 = ModernaStructure('file', fn, chain1)
            bp2 = ModernaStructure('file', fn, chain2)
            resi1 = [r for r in bp1][0]
            resi2 = [r for r in bp2][-1]
            hbonds = self.hb.calc_hbond_list(resi1, resi2)
            return hbonds

    def test_specific(self):
        # OK with just one hbond:
        # cHS_AC_Exemplar.pdb
        # tHS_UG_Exemplar.pdb
        #TODO: has to be argued:
        #  tWS_UU_Exemplar.pdb
        print(self.check_example('../../mod_isostericity/Westhof_new/tHS_AA_Exemplar.pdb'))
        
    def test_new_westhof_examples(self):
        wh_path = '../../mod_isostericity/Westhof_new/'
        for fn in os.listdir(wh_path):
            if not fn.endswith('.pdb'): continue
            hbonds = self.check_example(wh_path+fn)
            print(fn, hbonds)
          
            

if __name__ == '__main__':
    main()


