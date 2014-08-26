#!/usr/bin/env python
#
# test_check_pdb.py
#
# unit tests for the PDB cleanup module
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
from moderna.sequence.ModernaSequence import Sequence
from moderna.CheckPdb import PdbController
from moderna import load_model
import os
from test_data import *

class CheckPdbTests(TestCase):

    def setUp(self):
        self.st = ModernaStructure('file',NASTY_PDB)
        self.pc = PdbController(self.st)

    def test_ions_detection(self):
        self.assertEqual(len(self.pc.ions),9)

    def test_water_detection(self):
        self.assertEqual(len(self.pc.water),160)

    def test_unidentifiedRNA_detection(self):
        self.assertEqual(len(self.pc.unidentified_rna),1)

    def test_stars_detection(self):
        self.assertTrue(self.pc.stars)

    def test_O1P_detection(self):
        self.assertTrue(self.pc.OP)

    def test_clean_structure(self):
        self.pc.clean_structure()
        self.assertEqual(len(self.pc.ions),0)
        self.assertEqual(len(self.pc.ions),0)
        self.assertEqual(len(self.pc.unidentified_rna),1)
        self.assertEqual(len(self.pc.P_missing), 0)
        
        
    def test_clean_structure_in_struc(self):
        self.assertEqual(self.st.get_sequence(), Sequence('G_C_GGAU.UALCUCAGDDGGGAGAGCRCCAGABU#AAYAP?UGGAG7UC?UGUGTPCG"UCCACAGAAUUCGCACCA_._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.'))
        self.pc.clean_structure()
        self.assertEqual(self.st.get_sequence(), Sequence('GCGGAU.UALCUCAGDDGGGAGAGCRCCAGABU#AAYAP?UGGAG7UC?UGUGTPCG"UCCACAGAAUUCGCACCA'))
        
    def test_add_missing_p(self):
        """Identifies residues with P missing"""
        self.assertEqual(len(self.pc.P_missing), 1)
        self.assertEqual(self.pc.P_missing[0].identifier, '1')
        
    def test_atom_names(self):
        """Names of O1P, C1* etc should be fixed properly."""
        self.pc.clean_structure()
        at1 = self.st['2']["C1'"]
        self.assertEqual(at1.name, "C1'")
        self.assertEqual(at1.fullname, " C1'")
        at2 = self.st['4']["OP1"]
        self.assertEqual(at2.name, "OP1")
        self.assertEqual(at2.fullname, " OP1")
        at3 = self.st['2']["O5'"]
        self.assertEqual(at3.name, "O5'")
        self.assertEqual(at3.fullname, " O5'")
        
    def test_missing_phosphates(self):
        """residues with missing phosphates are fixed."""
        st2 = ModernaStructure('file',  MISSING_PHOSPHATES)
        pc = PdbController(st2)
        self.assertEqual(len(pc.P_missing), 2)
        self.assertEqual(st2.get_sequence(), Sequence('GCG_GAUUUALCUCAG'))
        pc.clean_structure()
        self.assertEqual(st2.get_sequence(), Sequence('GCG_GAUUUALCUCAG'))
        pc = PdbController(st2)
        self.assertEqual(len(pc.P_missing), 0)

    def test_missing_phosphates2(self):
        """Messy residues with missing phosphates are not delted."""
        st2 = ModernaStructure('file',  MISSING_PHOSPHATES2)
        pc = PdbController(st2)
        self.assertEqual(st2.get_sequence(), Sequence('C_C_G_A_C_C_U_U_C_G_G_C_C_A_C_C_U_G'))
        pc.clean_structure()
        self.assertEqual(st2.get_sequence(), Sequence('CCGACCUUCGGCCACC_UG'))

    def test_ligand_removal(self):
        """Should remove RNA ligand from 3FU2 when there is no ribose and phosphate group."""
        st2 = ModernaStructure('file',  PDB_WITH_LIGAND)
        pc2 = PdbController(st2)
        self.assertEqual(st2.get_sequence(), Sequence('AGAGGUUCUAG_._CACCCUCUAUAAAAAACUAA_x_._._._._._._._._._._._._.'))
        pc2.clean_structure()
        self.assertEqual(st2.get_sequence(), Sequence('AGAGGUUCUAG_CACCCUCUAUAAAAAACUAA'))
        
    def test_double_coordinates_identification(self):
        """Should remove RNA ligand when there is no ribose and phosphate group."""
        st3 = ModernaStructure('file',  PDB_WITH_DOUBLE_COORDINATES)
        pc3 = PdbController(st3)
        self.assertEqual(len(pc3.disordered), 3)
        pc3.clean_structure()
        self.assertEqual(len(pc3.disordered), 3)
        
    def test_AMP_detection_and_removal(self):
        """Should detect AMP in the structure and remove it while cleaning."""
        st4 = ModernaStructure('file',  PDB_WITH_AMP)
        self.assertEqual(st4.get_sequence(), Sequence('GGGAAGGGAAGAAACUGCGGCUUCGGCCGGCUUCCC_H'))
        pc4 = PdbController(st4)
        self.assertEqual(len(pc4.rna_ligand), 1)
        pc4.clean_structure()
        self.assertEqual(len(pc4.rna_ligand), 0)
        self.assertEqual(st4.get_sequence(), Sequence('GGGAAGGGAAGAAACUGCGGCUUCGGCCGGCUUCCC'))
        
    def test_phosphate_detection_and_removal(self):
        """Should detect pfosphate group (P, OP1, OP2, O5') in 3FU2 (resi 12) and remove it."""
        st5 = ModernaStructure('file',  PDB_WITH_PHOSPHATE)
        pc5 = PdbController(st5)
        self.assertEqual(len(pc5.phosphate), 1)
        pc5.clean_structure()
        self.assertEqual(len(pc5.phosphate), 0)
    
    def test_double_coordinates(self):
        """Disordered residues should be recognized only when they have double coordinates."""
        st6 = ModernaStructure('file',  PDB_WITHOUT_DOUBLE_COORD)
        pc6 = PdbController(st6)
        self.assertEqual(len(pc6.disordered), 0)
    
    
    
if __name__ == '__main__':
    main()
