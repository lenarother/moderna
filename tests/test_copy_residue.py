#!/usr/bin/env python
#
# test_copy_residue.py
#
# unit tests for residue duplication functionality
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
from moderna.ModernaResidue import ModernaResidue
from moderna.analyze.BaseRecognizer import BaseRecognizer
from moderna.util.Errors import ModernaStructureError

from test_data import *


class CopyResidueTests(TestCase):
    """
    Makes sure a ModernaResidue can be copied from
    one ModernaStructure to another (or itself),
    and that the resulting residue is really different.
    """
    def setUp(self):
        """Loads the A residue to start with."""
        self.struc_a = ModernaStructure('file',MINI_TEMPLATE)
        self.struc_b = ModernaStructure('file',MINI_TEMPLATE)
        self.struc_c = ModernaStructure()

    def test_copy(self):
        """The copied residue should appear at the new place."""
        copied = self.struc_a['5']
        self.assertEqual(copied.short_abbrev,'A')
        self.struc_c.add_residue(copied,'2')
        self.assertEqual(self.struc_c['2'].short_abbrev,'A')

    def test_non_identical(self):
        """Copy must produce a different instance."""
        self.struc_c.add_residue(self.struc_a['2'], '3')
        self.assertNotEqual(self.struc_c['3'],self.struc_a['2'])
        # KR: works only if there is no ModernaResidue.__cmp__

    def test_renumber(self):
        """New numbers should be applied automatically."""
        self.struc_c.add_residue(self.struc_a['2'],'3')
        self.assertEqual(self.struc_c['3'].identifier,'3')
        self.assertEqual(self.struc_a['2'].identifier,'2')

    def test_atom_number(self):
        """Copied residues should have the same number of atoms."""
        n_atoms_before = len(self.struc_a['5'])
        self.struc_c.add_residue(self.struc_a['5'],'30')
        self.assertEqual(len(self.struc_c['30']),n_atoms_before)
        
    def test_atom_non_identity(self):
        """Copied residues should have different atom instances."""
        self.struc_c.add_residue(self.struc_a['5'],'30')
        atoms_a = self.struc_a['5'].child_list
        atoms_c = self.struc_c['30'].child_list
        atoms_a.sort()
        atoms_c.sort()
        self.assertEqual(len(atoms_a),len(atoms_c))
        for i in range(len(atoms_a)):
            self.assertNotEqual(atoms_a[i], atoms_c[i])
        
    def test_atom_non_identity_in_structure(self):
        """Same residues from different structures should have different atom instances."""
        atoms_a = self.struc_a['5'].child_list
        atoms_b = self.struc_b['5'].child_list
        atoms_a.sort()
        atoms_b.sort()
        self.assertEqual(len(atoms_a),len(atoms_b))
        for i in range(len(atoms_a)):
            self.assertNotEqual(atoms_a[i], atoms_b[i])

    def test_modification(self):
        """Copying should also work for modified bases."""
        n_atoms_before = len(self.struc_a['10'])
        self.struc_b.add_residue(self.struc_a['10'],'30')
        self.assertEqual(len(self.struc_b['30']),n_atoms_before)
        
    def test_copy_exist(self):
        """ """
        self.struc_b.add_residue(self.struc_a['10'],'30')
        self.assertRaises(ModernaStructureError, self.struc_b.add_residue, self.struc_a['10'],'30')

        

if __name__ == '__main__':
    main()
  
