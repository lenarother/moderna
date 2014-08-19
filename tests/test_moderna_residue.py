#!/usr/bin/env python
#
# test_moderna_residue.py
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
from moderna.ModernaResidue import ModernaResidue
from moderna.analyze.BaseRecognizer import BaseRecognizer, BaseRecognitionError
from Bio.PDB import PDBParser
from moderna import load_model

from test_data import *

class ModernaResidueTests(TestCase):
    def setUp(self):
        """Loads the A residue to start with."""
        self.a=PDBParser().get_structure('test_struc',A_RESIDUE)[0].child_list[0].child_list[0]
        self.chain=PDBParser().get_structure('test_struc',MINI_TEMPLATE)[0].child_list[0]
        self.chain2=PDBParser().get_structure('test_struc',MINI_TEMPLATE)[0].child_list[0]
        self.unk=load_model(PDB_UNK)
    
    def tearDown(self):
        self.a = None
        self.chain = None
        self.chain2 = None

    def test_init(self):
        """Residues should be initializable."""
        a = ModernaResidue(self.a)
        self.assertEqual(a.identifier,'1')
        a = ModernaResidue(self.chain.child_list[0])
        self.assertEqual(a.identifier,'1')
        a = ModernaResidue(self.chain.child_list[-1])
        self.assertEqual(a.identifier,'15')
        
    def test_renumber(self):
        res = ModernaResidue(self.chain.child_list[2])
        res.change_number('64')
        self.assertEqual(res.identifier, '64')
        self.assertEqual(res.id, (' ',64, ' '))
        res.change_number('133E')
        self.assertEqual(res.identifier, '133E')
        self.assertEqual(res.id, (' ',133, 'E'))
        res.change_number('2')
        self.assertEqual(res.identifier, '2')
        self.assertEqual(res.id, (' ',2, ' '))

    def test_mutate(self):
        """Should allow all kinds of residue exchanges."""
        br = BaseRecognizer()
        resi = ModernaResidue(self.chain.child_list[2])
        abbrevs = ['A', 'Am', 'C', 'dC', 'dT', 'C','A', 'ac6A', 'm5U', 'Y', 'dA', 'Am', 'A']
        for name in abbrevs:
            resi.mutate(name)
            self.assertEqual(br.identify_resi(resi), name)
                         
    def test_mutate_unknown(self):
        """Should allow all kinds of residue exchanges."""
        br = BaseRecognizer()
        resi = ModernaResidue(self.chain.child_list[2])
        resi.mutate('X')
        self.assertRaises(BaseRecognitionError, br.identify_resi, resi)
        self.unk['39'].mutate('A')
        self.assertEqual(self.unk['39'].short_abbrev, 'A')
        self.unk['40'].mutate('Y')
        self.assertEqual(self.unk['40'].short_abbrev, 'P')

if __name__ == '__main__':
    main()
  
