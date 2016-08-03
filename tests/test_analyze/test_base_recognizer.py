#!/usr/bin/env python
"""
Unit Tests for the BaseRecognizer class
"""

from unittest import main, TestCase
from moderna.analyze.BaseRecognizer import BaseRecognizer
from Bio.PDB import PDBParser
from moderna.util.Errors import BaseRecognitionError
from moderna.Constants import PATH_TO_LIR_STRUCTURES
from test_data import *


class BaseRecognizerTests(TestCase):
    """
    Tests for the modification recognizer. Takes two sample
    PDB files and checks whether all modifications are found.
    """
    def setUp(self):
        self.br = BaseRecognizer()

    def tearDown(self):
        self.br = None

    def test_mini_template(self):
        """Should identify 15 bases including one modification."""
        names = ['G','C','G','G','A','U','U','U','A','m2G','C','U','C','A','G']
        struc = PDBParser().get_structure('test', MINI_TEMPLATE)
        chain = struc[0]['A']
        for resi, correct in zip(chain,names):
            base = self.br.identify_resi(resi)
            self.assertEqual(base,correct)
            
    def test_border_cases(self):
        """Recognize difficult residues by M.Skorupski."""
        path = TEST_DATA_PATH + 'nucleotides/border_cases/unknown%i_%s.pdb'
        EXAMPLES = [(1, 'C'), (2, 'C'), (3, 'C'), 
                    (4, 'U'), (5, 'C'), (6, 'C')]
        for num, base in EXAMPLES:
            fname = path % (num, base)
            struc = PDBParser().get_structure('test', fname)
            resi = struc[0].child_list[0].child_list[0]
            result = self.br.identify_resi(resi)
            self.assertEqual(result, base)

    def test_1ehz(self):
        """In the tRNA structure 14 modifications should be found."""
        # check the modifications in 1ehz
        ehz_modifications = dict([
            (10,'m2G'),(16,'D'),(17,'D'),(26,'m22G'),(32,'Cm'),
            (34,'Gm'),(37,'yW'),(39,'Y'),(40,'m5C'),(46,'m7G'),
            (49,'m5C'),(54,'m5U'),(55,'Y'),(58,"m1A")
            ])
        # merged with former speed test
        struc=PDBParser().get_structure('test',RNA_1EHZ)
        chain=struc[0]['A']
        checked = 0
        for resi in chain.child_list:
            base = self.br.identify_resi(resi)
            if ehz_modifications.has_key(resi.id[1]):
                self.assertEqual(base,ehz_modifications[resi.id[1]])
                checked += 1
            else:
                self.assertNotEqual(base,'')
        self.assertEqual(checked,14)


    def test_ms2i6A(self):
        """Test ms2i6A because this one is difficult."""
        struc=PDBParser().get_structure('test_struc',RNA_2OW8)
        chain=struc[0]['z']
        resi=chain.child_list[36]
        base = self.br.identify_resi(resi)
        self.assertEqual(base,'ms2i6A')

    def test_1qf6(self):
        """Test another difficult tRNA."""
        # check the modifications in 1qf6
        qf_modifications = dict([
            (16,'D'),(17,'D'),(20,'D'),(37,'m6t6A'),
            (46,'m7G'),(54,'m5U'),(55,'Y')
            ])
        struc=PDBParser().get_structure('test_struc', RNA_1QF6)
        chain=struc[0]['B']
        checked = 0
        for resi in chain.child_list:
            if qf_modifications.has_key(resi.id[1]):
                base = self.br.identify_resi(resi)
                self.assertEqual(base,qf_modifications[resi.id[1]])
                checked += 1
        self.assertEqual(checked,len(qf_modifications.keys()))
        
    def test_trna12H_37(self):
        """Difficult t6A example should work."""
        struc=PDBParser().get_structure('test_struc',PATH_TO_LIR_STRUCTURES+'trna12H.pdb')[0][' ']
        resi = struc.child_dict[('H_T6A',37,' ')]
        self.assertEqual(self.br.identify_resi(resi), 't6A')

    def test_pr0059Hc_32(self):
        """Difficult m2OH example should work."""
        struc=PDBParser().get_structure('test_struc',PATH_TO_LIR_STRUCTURES+'pr0059Hc.pdb')[0]['C']
        resi = struc.child_dict[('H_OMC',32,' ')]
        self.assertEqual(self.br.identify_resi(resi), 'Cm')
    
    def test_2ap0_9_a_8(self):
        """Difficult C example should work."""
        struc=PDBParser().get_structure('test_struc', TEST_DATA_PATH+'nucleotides/2ap0_9_A.pdb')[0]['A']
        resi = struc.child_dict[(' ',8,' ')]
        self.assertEqual(self.br.identify_resi(resi), 'C')
        
    def test_bromouridine(self):
        """Modifications with strange hetatoms should be caught."""
        struc=PDBParser().get_structure('test_struc', TEST_DATA_PATH+'nucleotides/5BrU.pdb')[0]['A']
        resi = struc.child_list[0]
        self.assertEqual(self.br.identify_resi(resi), '?U')
        
    def test_3jyv_7(self):
        """Difficult m2G example should work."""
        struc=PDBParser().get_structure('test_struc', TEST_DATA_PATH+'rna_structures/3jyv_7.pdb')[0]['7']
        resi = struc.child_dict[(' ',10,' ')]
        self.assertEqual(self.br.identify_resi(resi), 'm2G')
        resi = struc.child_dict[(' ',26,' ')]
        self.assertEqual(self.br.identify_resi(resi), 'm22G')
        resi = struc.child_dict[(' ',32,' ')]
        self.assertEqual(self.br.identify_resi(resi), 'Cm')
        resi = struc.child_dict[(' ',39,' ')]
        self.assertEqual(self.br.identify_resi(resi), 'Y')

    def test_atp_adp_amp(self):
        struc=PDBParser().get_structure('test_struc', TEST_DATA_PATH+'nucleotides/ATP_ADP_AMP.pdb')[0]['B']
        resis = struc.child_list
        self.assertEqual(self.br.identify_resi(resis[0]), 'ADP')
        self.assertEqual(self.br.identify_resi(resis[1]), 'GTP')
        self.assertEqual(self.br.identify_resi(resis[2]), 'AMP')
        self.assertEqual(self.br.identify_resi(resis[3]), 'GTP')
        self.assertEqual(self.br.identify_resi(resis[4]), 'dTMP')
        
    def test_amo(self):
        """recognize A with aminoacylated phosphate."""
        struc=PDBParser().get_structure('test_struc', TEST_DATA_PATH+'nucleotides/1C0A_AMO.pdb')[0]['A']
        resis = struc.child_list
        self.assertEqual(self.br.identify_resi(resis[0]), '?A_mod_phos')

    def test_so4_group(self):
        """Strange anorganic group should not work."""
        struc=PDBParser().get_structure('test_struc',SULFATE)[0]['A']
        resi = struc.child_dict[('H_SO4',1481,' ')]
        self.assertRaises(BaseRecognitionError,self.br.identify_resi,resi)
        
    def test_misc_ligands(self):
        """Large set of exemplars is recognized."""
        expected = [
            '3pADP', 'GMP', 'm66A', 'ms2i6Aiso', 'A', 'd8fG', 'preQ1tRNA', '23pC', 'D', 'arabinoseU', 
            'd8oG', '2ofluoro-m5U', '3meo5mC', 'Arp', 'dmh5U', 'QtRNA', 'm5Ueth', 'dG', 'm1A', '2oNAcC', 
            'm4Cm', '3meoG', 'phosphonoG', 'm6t6A', '2ofluoro-8oG', 'd5propnU', 'm1A', 'd5mpC', 'm1G', 'm7G', 
            'mnm5U', 'o6U', 'de3T', 'm22G', 'm7G', "?Tm_5'amino", 'disoG', '?U_mod_phos', 'm5U', 'mcm5s2U', 
            '2oNC', 'dA', 'arabinoseA', '7fluorobenzylG', 'm5C', 'VA', 'd3pT', 'mnm5s2U', 'dhOroP', 'm5Um', 
            'd5mpA', '5fluoroC', 's2U', 'dm5C', '3meoA', 'Am', '?C_mod_phos', 'd5mpA', 'd3pT', '5fluoroU', 
            'dG', 'd5propU', 'N2A', 'alpha-dA', 'dC', '?G_mod_phos', '2ofluoro-m7G', 'd5nitroU', 'dN2A', 'Cm', 
            'dethC', '2mcarbT', '23pA', 'Gm', 'D', 'yW', '5nitroC', 'Um', 's4U', None, 
            '?U_mod_diphos', 'm5Um', 'arabinoseC', 'm2A', 'm2G', 't6A', 'd35pG', 'A', 'dm8G', 'm5U', 
            'arabinoseG', None, '35pG', '?A_mod_phos', '?U', 'dm6G', 'tfA', '?U_mod_diphos', '?A_mod_phos', 'dmo5U', 
            None, 'dN72G', '3meo5mU', 'tfT', 'Y', '2ofluoro-8oG', "?Tm_5'amino", 'phosphonoG', 'd5propU'
        ]
        correct = 0
        chain=PDBParser().get_structure('test_struc', TEST_DATA_PATH+'nucleotides/misc_nucleotides.pdb')[0]['A']
        for resi, exp in zip(chain.child_list, expected):
            try:
                result = self.br.identify_resi(resi)
            except BaseRecognitionError:
                result = None
            if result == exp:
                correct += 1
            else:
                print(resi, result, exp)
        self.assertEqual(correct, len(chain.child_list))

    def test_dna(self):
        """DNA should also be recognized correctly."""
        struc = PDBParser().get_structure('test_struc',DNA_WITH_MISMATCH)
        chain = struc[0]['E']
        seq = [self.br.identify_resi(resi) for resi in chain.child_list]
        self.assertEqual(seq, ['dA','dG','dC','dT','dG','dC','dC','dA','dG',\
                               'dG','dC','dA','dC','dC','dA','dG','dT','dG'])

    def test_fail_unrecognizable(self):
        """Bad residues should raise an exception."""
        struc=PDBParser().get_structure('test_struc',UNRECOGNIZABLE_RESIDUE)
        resi = struc[0]['A'].child_list[0]
        self.assertRaises(BaseRecognitionError, self.br.identify_resi, resi)
        
    def test_too_short_glyc(self):
        """dA with too short glycosidic bond should be recognized."""
        struc=PDBParser().get_structure('test_struc',A_TOO_SHORT)
        resi = struc[0]['A'].child_list[0]
        self.assertEqual(self.br.identify_resi(resi), 'dA')
        
    def test_richardson(self):
        """Should get sequence of a Richardson RNADB2005 file correctly."""
        struc=PDBParser().get_structure('test_struc',RICHARDSON_STRUC)
        chain=struc[0]['A']
        SEQ = "UUAUAUAUAUAUAA"
        for resi, expected in zip(chain.child_list, SEQ):
            self.assertEqual(self.br.identify_resi(resi), expected)

    def test_protein(self):
        """Amino acids should be recognized as well."""
        # check the modifications in 1ehz
        amino = dict([
            (1,'ALA'), (2, 'GLN'), (3, 'THR'),  (4, 'VAL') , (5, 'PRO'), 
            (6, 'TYR'), (7, 'GLY'), (8, 'ILE'), (10, 'LEU'), (12, 'LYS'), 
            (14, 'ASP'),(21, 'PHE') , (38, 'SER'), (39, 'HIS'), (43, 'ASN'), 
            (119, 'MET'), (113,'TRP'), (54, 'GLU') # arg and trp missing
            ])
        struc=PDBParser().get_structure('test_struc',PROTEIN_STRUCTURE)
        chain=struc[0]['E']
        checked = 0
        for resi in chain.child_list:
            if amino.has_key(resi.id[1]):
                aa = self.br.identify_resi(resi)
                self.assertEqual(aa,amino[resi.id[1]])
                checked += 1
        self.assertEqual(checked,18)
        
    def test_cap5(self):
        """Caps in mRNA are recognized."""
        struc=PDBParser().get_structure('test_struc', CAP_EXAMPLE)
        resi = struc[0]['B'].child_list[0]
        self.assertEqual(self.br.identify_resi(resi), 'm7Gpp_cap')



if __name__ == '__main__':
    main()
    
