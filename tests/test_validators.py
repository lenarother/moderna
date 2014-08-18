#!/usr/bin/env python
#
# test_validators.py
#
# unit tests for parameter validators.
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
from moderna.Errors import ParameterError
from moderna.ModernaStructure import ModernaStructure
from moderna.Template import Template
from moderna.RNAModel import RnaModel
from moderna.ModernaFragment import ModernaFragment53
from moderna.ModernaResidue import ModernaResidue
from moderna.sequence.RNAAlignment import read_alignment
from moderna.sequence.ModernaSequence import Sequence
from Bio.PDB import PDBParser
from test_data import TEST_DATA_PATH, A_RESIDUE, MINI_TEMPLATE, MINI_ALIGNMENT_FILE
from moderna.validators import validate_structure, validate_template, validate_model, \
    validate_alignment, validate_seq, validate_resnum, \
    validate_resi, validate_fragment, \
    validate_resi_list, validate_filename, validate_path, \
    validate_secstruc, validate_alphabet, validate_alphabet_list


class StrucValidatorTests(TestCase):
    """
    Tests for parameter validation functions with structures.
    """
    def setUp(self):
        self.s = ModernaStructure('file', MINI_TEMPLATE)
        self.t = Template(MINI_TEMPLATE, 'file')
        self.m = RnaModel(None, None, 'A', 'file', MINI_TEMPLATE)
        self.f = ModernaFragment53(self.s, anchor5=self.m['1'], anchor3=self.m['14'])
        self.struc=PDBParser().get_structure('test_struc', MINI_TEMPLATE)
        
    def test_validate_structure(self):
        """Checks for ModernaStructure objects."""
        self.assertEqual(validate_structure(self.s), self.s)
        self.assertEqual(validate_structure(self.t), self.t)
        self.assertEqual(validate_structure(self.m), self.m)
        self.assertRaises(ParameterError, validate_model, self.f)
        self.assertRaises(ParameterError, validate_structure, self.struc)
        
    def test_validate_template(self):
        self.assertEqual(validate_template(self.t), self.t)
        self.assertRaises(ParameterError, validate_template, self.s)
        self.assertRaises(ParameterError, validate_template, self.m)
        self.assertRaises(ParameterError, validate_model, self.f)
        self.assertRaises(ParameterError, validate_template, self.struc)

    def test_validate_model(self):
        self.assertEqual(validate_model(self.m), self.m)
        self.assertRaises(ParameterError, validate_model, self.s)
        self.assertRaises(ParameterError, validate_model, self.t)
        self.assertRaises(ParameterError, validate_model, self.f)
        self.assertRaises(ParameterError, validate_model, self.struc)
        
    def test_validate_fragment(self):
        self.assertEqual(validate_fragment(self.f), self.f)
        self.assertRaises(ParameterError, validate_fragment, self.s)
        self.assertRaises(ParameterError, validate_fragment, self.t)
        self.assertRaises(ParameterError, validate_fragment, self.m)
        self.assertRaises(ParameterError, validate_fragment, self.struc)


class ValidatorTests(TestCase):
    """
    Tests for parameter validation functions.
    """
    def test_validate_alignment(self):
        ali_str = "> target\n--AAAA\n> template\nGGGG--\n"
        ali = read_alignment(ali_str)
        exp_ali = read_alignment(MINI_ALIGNMENT_FILE)
        self.assertEqual(str(validate_alignment(ali)), str(ali))
        self.assertEqual(str(validate_alignment(ali_str)), str(ali))
        self.assertEqual(str(validate_alignment(MINI_ALIGNMENT_FILE)), str(exp_ali))
        self.assertRaises(ParameterError, validate_alignment, 112345)
        
    def test_validate_seq(self):
        seq_str = "AGCU:PY_7"
        seq = Sequence(seq_str)
        self.assertEqual(validate_seq(seq), seq)
        self.assertEqual(validate_seq(seq_str), seq)
        self.assertEqual(validate_seq(''), Sequence(''))
        self.assertRaises(ParameterError, validate_seq, 12345)

    def test_validate_secstruc(self):
        self.assertEqual(validate_secstruc('((...))'), '((...))')
        self.assertEqual(validate_secstruc('...((('), '...(((')
        self.assertEqual(validate_secstruc(''), '')
        self.assertRaises(ParameterError, validate_secstruc, '[[]]')
        self.assertRaises(ParameterError, validate_secstruc, '{{}}')
        self.assertRaises(ParameterError, validate_secstruc, '<<>>')

    def test_validate_resnum(self):
        self.assertEqual(validate_resnum('5'), '5')
        self.assertEqual(validate_resnum(123), '123')
        self.assertEqual(validate_resnum('27G'), '27G')
        self.assertRaises(ParameterError, validate_resnum, 'A')
        self.assertRaises(ParameterError, validate_resnum, 'A1')
        self.assertRaises(ParameterError, validate_resnum, ' 1 ')
        self.assertRaises(ParameterError, validate_resnum, '')

    def test_validate_alphabet(self):
        self.assertEqual(validate_alphabet('m1A'), 'm1A')
        self.assertEqual(validate_alphabet('7'), 'm7G')
        self.assertRaises(ParameterError, validate_alphabet, '999U')

    def test_validate_alphabet_list(self):
        self.assertEqual(validate_alphabet_list(['m1A', 'G']), ['m1A', 'G'])
        self.assertEqual(validate_alphabet_list('AG'), ['A', 'G'])
        self.assertRaises(ParameterError, validate_alphabet, 123)
        self.assertRaises(ParameterError, validate_alphabet, ['A', '_'])

    def test_validate_resi(self):
        m = ModernaStructure('file', MINI_TEMPLATE)
        self.assertEqual(validate_resi(m['5']), m['5'])
        self.assertRaises(ParameterError, validate_resi, m)
        self.assertRaises(ParameterError, validate_resi, 'not a residue')

    def test_validate_resi_list(self):
        m = ModernaStructure('file', MINI_TEMPLATE)
        expected = list(m)
        self.assertEqual(validate_resi_list(m), expected)
        self.assertEqual(validate_resi_list([m['2'], m['3']]), expected[1:3])
        self.assertRaises(ParameterError, validate_resi_list, ['not a residue'])
        self.assertRaises(ParameterError, validate_resi_list, 123)
        
    def test_validate_filename(self):
        self.assertEqual(validate_filename(TEST_DATA_PATH+'nonexist.bla'), TEST_DATA_PATH+'nonexist.bla')
        self.assertEqual(validate_filename(A_RESIDUE, True), A_RESIDUE)
        self.assertRaises(ParameterError, validate_filename, TEST_DATA_PATH+'nonexist.bla', True)
        
    def test_validate_path(self):
        self.assertEqual(validate_path(''), '')
        self.assertEqual(validate_path(TEST_DATA_PATH), TEST_DATA_PATH)
        self.assertRaises(ParameterError, validate_path, '/@#%$@#%\t@#forsurethere/isnosuch/dir^%()^\n')
    
if __name__ == '__main__':
    main()
