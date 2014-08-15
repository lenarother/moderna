
from moderna.renumber_chain import renumber_chain, renumber_all_residues, insert_gap
from moderna.ModernaStructure import ModernaStructure
from moderna.Errors import RenumeratorError
from unittest import TestCase, main
from test_data import RNA_1EHZ, MINI_TEMPLATE


class RenumberChainTests(TestCase):
    """
    Tests for the abstract renumerator class
    """
    def setUp(self):
        self.rna = ModernaStructure('file', MINI_TEMPLATE)
    
    def test_renumber_all_residues(self):
        """Renumber chain with a list of numbers"""
        numbers = map(str, range(1, 16))
        renumber_all_residues(self.rna, numbers)
        self.assertEqual(len(self.rna), 15)
        self.assertEqual(self.rna['1'].long_abbrev, 'G')
        self.assertEqual(self.rna['15'].long_abbrev, 'G')
        
    def test_renumber_all_residues_error(self):
        """Too short list of numbers raises exception"""
        numbers = map(str, range(1, 15))
        self.assertRaises(RenumeratorError, renumber_all_residues, self.rna, numbers)
        
    def test_renumber_chain(self):
        """Renumber entire chain."""
        m = ModernaStructure('file', RNA_1EHZ)
        renumber_chain(m, '100')
        self.assertEqual(m['100'].long_abbrev, 'G')
        self.assertEqual(m['107'].long_abbrev, 'U')
        self.assertEqual(m['109'].long_abbrev, 'm2G')
        self.assertEqual(m['115'].long_abbrev, 'D')
        self.assertEqual(m['175'].long_abbrev, 'A')
        
    def test_renumber_chain_default(self):
        """Default starting number is 1"""
        renumber_chain(self.rna)
        self.assertEqual(self.rna['1'].long_abbrev, 'G')
        self.assertEqual(self.rna['15'].long_abbrev, 'G')
    
    def test_renumber_chain_error(self):
        """Non-integer numbers give an error"""
        self.assertRaises(RenumeratorError, renumber_chain, self.rna, '77A')
        
    def test_insert_gap(self):
        """insert gap that shifts numbers"""
        insert_gap(self.rna, '5', 5)
        numbers = [resi.identifier for resi in self.rna]
        self.assertEqual(numbers, ['1', '2', '3', '4', '5', '11', '12',\
             '13', '14', '15', '16', '17', '18', '19', '20'])

    def test_insert_gap_integer(self):
        """insert gap accepts integer input"""
        insert_gap(self.rna, 5, 5)

    def test_insert_gap_none(self):
        """insert gap accepts none value aka start"""
        insert_gap(self.rna, None, 5)
        numbers = [resi.identifier for resi in self.rna]
        self.assertEqual(numbers, ['6', '7', '8', '9', '10', '11', '12',\
             '13', '14', '15', '16', '17', '18', '19', '20'])

    def test_insert_gap_twice(self):
        """double gap inserts just enough numbers"""
        insert_gap(self.rna, '5', 5)
        insert_gap(self.rna, '5', 7)
        numbers = [resi.identifier for resi in self.rna]
        self.assertEqual(numbers, ['1', '2', '3', '4', '5',\
             '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'])
        


if __name__ == '__main__':
    main()
    
