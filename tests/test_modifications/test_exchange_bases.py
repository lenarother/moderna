#!/usr/bin/env python
#
# test_exchange_bases.py
#
# unit tests for base exchange functionality
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
from moderna.RNAResidue import RNAResidue
from moderna.analyze.BaseRecognizer import BaseRecognizer
from moderna.modifications import exchange_base, make_backbone_only_residue, modify_residue
from moderna.util.Errors import BaseRecognitionError
from test_data import *
from Bio.PDB.Vector import calc_dihedral
from Bio.PDB.PDBParser import PDBParser


class ExchangeBaseTests(TestCase):
    """
    Makes sure the four standard bases A,G,C,U
    can be replaced in ModernaResidues by each
    other.
    """
    def setUp(self):
        """Loads the A residue to start with."""
        self.struc = ModernaStructure('file',A_RESIDUE)
        self.adenosine = self.struc['1']
        struc = ModernaStructure('file',C_RESIDUE, '0')
        self.cytosine = struc['494']

    def test_init(self):
        """The Adenosine test file must be really an A."""
        self.assertEqual(self.adenosine.long_abbrev,'A')
        self.assertEqual(self.adenosine.short_abbrev,'A')
        recon = BaseRecognizer().identify_resi(self.adenosine)
        self.assertEqual(recon,'A')

    def test_substitute(self):
        """Exchange must produce G,C, and U."""
        exchange_base(self.adenosine, 'G')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'G')
        exchange_base(self.adenosine, 'C')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'C')
        exchange_base(self.adenosine, 'U')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'U')

    def test_resubstitute(self):
        """Substituting A by U and U by A should give A again."""
        exchange_base(self.adenosine, 'G')
        exchange_base(self.adenosine, 'A')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'A')
        #
        exchange_base(self.adenosine, 'C')
        exchange_base(self.adenosine, 'A')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'A')
        #
        exchange_base(self.adenosine, 'U')
        exchange_base(self.adenosine, 'A')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'A')

    def test_same_dihedral(self):
        """The dihedral of the glycosidic bond before and after
        substitution should be the same."""
        # error in dihedral angle is 1.0 degrees.
        # before (A)
        a1 = self.adenosine["C2'"].get_vector()
        a2 = self.adenosine["C1'"].get_vector()
        a3 = self.adenosine["N9"].get_vector()
        a4 = self.adenosine["C4"].get_vector()
        a_dihedral = calc_dihedral(a1,a2,a3,a4)
        # for G
        exchange_base(self.adenosine, 'G')
        g = self.adenosine                                     
        g1 = g["C2'"].get_vector()
        g2 = g["C1'"].get_vector()
        g3 = g["N9"].get_vector()
        g4 = g["C4"].get_vector()
        g_dihedral = calc_dihedral(g1,g2,g3,g4)
        self.assertAlmostEqual(g_dihedral,a_dihedral,2)
        # for U
        exchange_base(self.adenosine, 'U')
        u = self.adenosine
        u1 = u["C2'"].get_vector()
        u2 = u["C1'"].get_vector()
        u3 = u["N1"].get_vector()
        u4 = u["C2"].get_vector()
        u_dihedral = calc_dihedral(u1,u2,u3,u4)
        self.assertAlmostEqual(u_dihedral,a_dihedral,2)

class ModifyResidueTests(TestCase):

    def setUp(self):
        """Loads the A residue to start with."""
        self.struc = ModernaStructure('file',A_RESIDUE)
        self.adenosine = self.struc['1']
        struc = ModernaStructure('file',C_RESIDUE, '0')
        self.cytosine = struc['494']
        self.chain = PDBParser().get_structure('test_struc',MINI_TEMPLATE)[0].child_list[0]
        self.unk = ModernaStructure('file', PDB_UNK)

    def test_make_backbone_only_residue(self):
        """Should leave only backbone and ribose"""
        make_backbone_only_residue(self.adenosine)
        atomnames = "P,O5',C5',C4',C3',O3',OP1,OP2,O2',C2',C1',O4',N9".split(',')
        self.assertEqual(len(self.adenosine.child_list), len(atomnames))
        for atom in self.adenosine.child_list:
            self.assertTrue(atom.name in atomnames)
        
        # second example: cytosine
        make_backbone_only_residue(self.cytosine)
        atomnames = "P,O5',C5',C4',C3',O3',OP1,OP2,O2',C2',C1',O4',N1".split(',')
        self.assertEqual(len(self.cytosine.child_list), len(atomnames))
        for atom in self.cytosine.child_list:
            self.assertTrue(atom.name in atomnames)
            
    def test_backbone_getitem(self):
        """Backbone residues should have a N*"""
        make_backbone_only_residue(self.adenosine)
        self.assertTrue(self.adenosine['N*'])
        make_backbone_only_residue(self.cytosine)
        self.assertTrue(self.cytosine['N*'])
        

    def test_mutate(self):
        """Should allow all kinds of residue exchanges."""
        br = BaseRecognizer()
        resi = RNAResidue(self.chain.child_list[2])
        abbrevs = ['A', 'Am', 'C', 'dC', 'dT', 'C','A', 'ac6A', 'm5U', 'Y', 'dA', 'Am', 'A']
        for name in abbrevs:
            modify_residue(resi, name)
            self.assertEqual(br.identify_resi(resi), name)

    def test_mutate_unknown(self):
        """Should allow all kinds of residue exchanges."""
        br = BaseRecognizer()
        resi = RNAResidue(self.chain.child_list[2])
        modify_residue(resi, 'X')
        self.assertRaises(BaseRecognitionError, br.identify_resi, resi)
        modify_residue(self.unk['39'], 'A')
        self.assertEqual(self.unk['39'].short_abbrev, 'A')
        modify_residue(self.unk['40'], 'Y')
        self.assertEqual(self.unk['40'].short_abbrev, 'P')

if __name__ == '__main__':
    main()
  
