
from moderna.Renumerator import Renumerator, NumberRenumerator, LetterRenumerator
from unittest import TestCase, main
from moderna.Errors import RenumeratorError


class RenumeratorTests(TestCase):
    """
    Tests for the abstract renumerator class
    """
    def test_divide_identifiers(self):
        """Splits a list of identifiers."""
        r = Renumerator([None, None, None, '2', '3'])
        self.assertEqual(str(r.divide_identifiers()), "[[None, None, None], ['2', '3']]")
        r = Renumerator(['5', '6', None, None, None, '7', '8'])
        self.assertEqual(str(r.divide_identifiers()), "[['5', '6'], [None, None, None], ['7', '8']]")
        r = Renumerator(['1', '2', None, None, None, '3', '4', '5', None, None, None, '6', '7'])
        self.assertEqual(str(r.divide_identifiers()),"[['1', '2'], [None, None, None], ['3', '4', '5'], [None, None, None], ['6', '7']]")
    
    

class NumberRenumeratorTests(TestCase):

    def test_number_generator(self):
        """Increase a number stepwise"""
        r = NumberRenumerator(['5', '6', None, None, None, '7', '8'])
        gen = r.number_generator('42')
        self.assertEqual(gen.next(), '43')
        self.assertEqual(gen.next(), '44')
        self.assertEqual(gen.next(), '45')

    def test_numbers_from_scratch(self):
        r = NumberRenumerator([None, None, None])
        result = r.get_identifiers()
        self.assertEqual(result, ['1', '2', '3'])
    
    def test_5_prime(self):
        r = NumberRenumerator([None, None, None, '2', '3'])
        result = r.get_identifiers()
        self.assertEqual(result, ['1', '2', '3', '4', '5'])

    def test_5_prime_long(self):
        r = NumberRenumerator([None]*50 + ['2', '3'])
        result = r.get_identifiers()
        self.assertEqual(result[-1], '52')
        r = NumberRenumerator([None]*50 + ['200', '201'])
        result = r.get_identifiers()
        self.assertEqual(result[0], '150')

    def test_3_prime(self):
        r = NumberRenumerator(['5', '6', None, None, None])
        result = r.get_identifiers()
        self.assertEqual(result, ['5', '6', '7', '8', '9'])

    def test_53_insert(self):
        r = NumberRenumerator(['5', '6', None, None, None, '7', '8'])
        result = r.get_identifiers()
        self.assertEqual(result, ['5', '6', '7', '8', '9', '10', '11'])

    def test_twin_segments(self):
        r = NumberRenumerator(['1', '2', None, None, None, '3', '4', '5', None, None, None, '6', '7'])
        result = r.get_identifiers()
        self.assertEqual(result, ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13'])
        

class LetterRenumeratorTests(TestCase):
    """
    Tests for the renumerator class
    """
    def test_letter_generator(self):
        """Append to a number A, B, C, ..."""
        r = LetterRenumerator(['5', '6', None, None, None, '7', '8'])
        gen = r.letter_generator('42')
        self.assertEqual(gen.next(), '42A')
        self.assertEqual(gen.next(), '42B')
        self.assertEqual(gen.next(), '42C')
    
    def test_numbers_from_scratch(self):
        r = LetterRenumerator([None, None, None])
        result = r.get_identifiers()
        self.assertEqual(result, ['1A', '1B', '1C'])
    
    def test_5_prime(self):
        r = LetterRenumerator([None, None, None, '2', '3'])
        result = r.get_identifiers()
        self.assertEqual(result, ['1A', '1B', '1C', '2', '3'])
        
    def test_5_prime_long(self):
        r = LetterRenumerator([None]*50 + ['2', '3'])
        self.assertRaises(RenumeratorError, r.get_identifiers)
        r = LetterRenumerator([None]*50 + ['200', '201'])
        result = r.get_identifiers()
        self.assertEqual(result[0], '150')

    def test_3_prime(self):
        r = LetterRenumerator(['5', '6', None, None, None])
        result = r.get_identifiers()
        self.assertEqual(result, ['5', '6', '6A', '6B', '6C'])
        
    def test_53_insert(self):
        r = LetterRenumerator(['5', '6', None, None, None, '7', '8'])
        result = r.get_identifiers()
        self.assertEqual(result, ['5', '6', '6A', '6B', '6C', '7', '8'])
        
    def test_twin_segments(self):
        r = LetterRenumerator(['1', '2', None, None, None, '3', '4', '5', None, None, None, '6', '7'])
        result = r.get_identifiers()
        self.assertEqual(result, ['1', '2', '2A', '2B', '2C', '3', '4', '5', '5A', '5B', '5C', '6', '7'])
   
    def test_prepare_identifiers(self):
        """Should generate a list of identifiers."""
        r = LetterRenumerator([])
        # letters
        result = r.prepare_identifiers(['10'], [None]*5, ['20'])
        self.assertEqual(result, ['10A', '10B', '10C', '10D', '10E'])
        # numbers
        result = r.prepare_identifiers(['10'], [None]*50, ['100'])
        self.assertEqual(result[:5], ['11', '12', '13', '14', '15'])
        # Exception when there is not enough room
        self.assertRaises(RenumeratorError, r.prepare_identifiers, ['10'], [None]*50, ['20'])


if __name__ == '__main__':
    main()
    
