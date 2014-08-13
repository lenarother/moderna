#!/usr/bin/env python
#
# test_write_pdb.py
#
# tests for ModernaStructure PDB write feature.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

from unittest import main, TestCase
from moderna.ModernaStructure import ModernaStructure
from Bio.PDB import PDBParser
from test_data import *

TEST_OUTPUT = 'test_data/test_output.ent'

class WritePDBTests(TestCase):

    def get_resnames(self,fn,chain='A'):
        """Returns a list of residue names from a PDB file."""
        result = []
        struc=PDBParser().get_structure('test_struc',fn)
        chain=struc[0][chain]
        for resi in chain.child_list:
            result.append(resi.resname)
        return result

    def setUp(self):
        self.s = ModernaStructure('file',MINI_TEMPLATE)

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
        self.s['5'].exchange_base('C')
        self.s.write_pdb_file(TEST_OUTPUT)
        rn = self.get_resnames(TEST_OUTPUT)
        self.assertEqual(rn[4], '  C')
        
    def test_res_names_add_modif(self):
        """Names should be consistent after adding modifications."""
        self.s['5'].add_modification('m7G')
        self.s.write_pdb_file(TEST_OUTPUT)
        rn = self.get_resnames(TEST_OUTPUT)
        self.assertEqual(rn[4], '7MG')
        
    def test_res_names_remove_modif(self):
        """Names should be consistent after adding modifications."""
        self.s['10'].remove_modification()
        self.s.write_pdb_file(TEST_OUTPUT)
        rn = self.get_resnames(TEST_OUTPUT)
        self.assertEqual(rn[4], '  A')

    def test_atom_number(self):
        pass

if __name__ == '__main__':
    main()
    
        
