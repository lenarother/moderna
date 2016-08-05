#!/usr/bin/env python
"""
tests for ModernaStructure PDB write feature.
"""

from unittest import main, TestCase
from moderna.ModernaStructure import ModernaStructure
from moderna.modifications import exchange_base, add_modification, remove_modification
from Bio.PDB import PDBParser
from moderna.tests.test_data import *
import os

TEST_OUTPUT = 'test_output.ent'


class WritePDBTests(TestCase):

    def get_resnames(self,fn,chain='A'):
        """Returns a list of residue names from a PDB file."""
        result = []
        struc=PDBParser().get_structure('test_struc',fn)
        chain=struc[0][chain]
        for resi in chain.child_list:
            result.append(resi.resname)
        return result

    def remove_outfile(self):
        if os.path.exists(TEST_OUTPUT):
            os.remove(TEST_OUTPUT)

    def setUp(self):
        self.remove_outfile()
        self.s = ModernaStructure('file',MINI_TEMPLATE)

    def tearDown(self):
        self.remove_outfile()

    def test_res_names(self):
        """Names of written residues should be standard."""
        self.s.write_pdb_file(TEST_OUTPUT)
        rn = self.get_resnames(TEST_OUTPUT)
        self.assertEqual(rn, ['  G', '  C', '  G', '  G', '  A', '  U', '  U', '  U', '  A', '2MG', '  C', '  U', '  C', '  A', '  G'])

    def test_res_names_mod(self):
        """Names of modifications should be standard."""
        # I'll check file with amber names manualy. 
        # There are some conflicts and inconsistents.
        pass
        
    def test_res_names_exchanged(self):
        """Names should be consistent after base exchanges."""
        exchange_base(self.s['5'], 'C')
        self.s.write_pdb_file(TEST_OUTPUT)
        rn = self.get_resnames(TEST_OUTPUT)
        self.assertEqual(rn[4], '  C')
        
    def test_res_names_add_modif(self):
        """Names should be consistent after adding modifications."""
        add_modification(self.s['5'], 'm7G')
        self.s.write_pdb_file(TEST_OUTPUT)
        rn = self.get_resnames(TEST_OUTPUT)
        self.assertEqual(rn[4], '7MG')
        
    def test_res_names_remove_modif(self):
        """Names should be consistent after adding modifications."""
        remove_modification(self.s['10'])
        self.s.write_pdb_file(TEST_OUTPUT)
        rn = self.get_resnames(TEST_OUTPUT)
        self.assertEqual(rn[4], '  A')

    def test_atom_number(self):
        pass

if __name__ == '__main__':
    main()
    
        
