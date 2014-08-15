#!/usr/bin/env python
#
# Renumerator.py
#
# Generates identifiers for residues from ModernaFragment objects.
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


from Errors import RenumeratorError

class Renumerator:
    """
    Class for generating number for residues renumeration.
    Useful in fragments of all kinds
    """
    def __init__ (self, identifiers):
        self.identifiers = identifiers

    def divide_identifiers(self):
        """
        divides the list of input identifiers
        e.g.
        [None, None, None, '2', '3'] ---> [[None, None, None], ['2', '3']]
        ['5', '6', None, None, None, '7', '8'] ---> [['5', '6'], [None, None, None], ['7', '8']]
        ['1', '2', None, None, None, '3', '4', '5', None, None, None, '6', '7'] ---> [['1', '2'], [None, None, None], ['3', '4', '5'], [None. None, None], ['6', '7']]
        """
        divided = []
        num = self.identifiers[0]
        temp = []
        for x in self.identifiers:
            if (x and num) or (not x and not num):
                temp.append(x)
            else:
                divided.append(temp)
                temp = [x]
            num = x
        if temp: 
            divided.append(temp)
        return divided
        
    def get_generator(self, first, before, middle, after):
        '''
        abstract method that generates identifiers.
        '''
        raise RenumeratorError('Abstract Method - use a Renumerator subclass!')

    def prepare_identifiers(self, before, middle, after):
        """
        before - list with identifiers before query, or None
        middle - list with None elements (query)
        after - list with identifiers after query or None
        """
        fn = self.find_first_number(before, middle, after)
        generator = self.get_generator(fn, before, middle, after)
        identifiers = [generator.next() for x in middle]
        return identifiers

    def get_sections(self, divided, x):
        before, after = None, None
        if x > 0:
            before = divided[x-1]
        if x < len(divided)-1: 
            after = divided[x+1]
        return before, divided[x], after
        
    def get_identifiers(self):
        """
        Returns a complete list of new identifiers for the given list of identifiers.
        """
        identifiers = []
        divided = self.divide_identifiers()
        for x in range(len(divided)):
            before, middle, after = self.get_sections(divided, x)
            if middle[0] == None:
                identifiers += self.prepare_identifiers(before, middle, after)
            else: 
                identifiers += middle
        return identifiers


class NumberRenumerator(Renumerator):
    '''
    Assigns new numbers as integer numbers.
    Changes existing numbers as well if necessary.
    
    e.g.
      [None, None, None] ---> ['1A', '1B', '1C']
      [None, None, None, '2', '3'] ---> ['1A', '1B', '1C', '2', '3']
      ['5', '6', None, None, None, '7', '8'] ---> ['5', '6', '6A', '6B', '6C', '7', '8']
      ['1', '2', None, None, None, '3', '4', '5', None, None, None, '6', '7'] ---> ['1', '2', '2A', '2B', '2C', '3', '4', '5', '5A', '5B', '5C', '6', '7']
     '''
    def number_generator(self, first_num):
        try: 
            first_num = int(first_num)
        except ValueError: 
            RenumeratorError('Number can not be generated, first number cannot be converted into int')
        while 1:
            first_num += 1
            yield str(first_num)

    def find_first_number(self, before, middle, after):
        """
        Returns starting number.
        """
        if not before and not after:
            return '0'
        elif before:
            return before[-1]
        else: 
            start = int(after[0]) - len(middle) -1 
            if start < 0:
                return '0'
            return start
            
    def get_generator(self, first, before, middle, after):
        return self.number_generator(first)
 
    def get_identifiers(self):
        """
        Returns a complete list of new identifiers for the given list of identifiers.
        """
        identifiers = []
        divided = self.divide_identifiers()
        
        for x in range(len(divided)):
            before, middle, after = self.get_sections(divided, x)
            if len(identifiers) > 0:
                lastid = identifiers[-1]
            else:
                lastid = '0'
            if before and before[0] != None and int(before[-1]) < int(lastid):
                before = [lastid]
            if middle[0] == None:
                identifiers += self.prepare_identifiers(before, middle, after)
            elif int(lastid) >= int(middle[0]):
                identifiers += self.prepare_identifiers([lastid], middle, None)
            else:
                identifiers += middle
        return identifiers        
        
class LetterRenumerator(Renumerator):  
    '''
    Assigns new numbers as insertion codes (3A, 3B,..).
    Never changes existing numbers.
    '''
    def letter_generator(self, first_num): 
        CHARACTERS = 'ABCDEFGHIJKLMNOPQRSTUWXYZ'
        for x in CHARACTERS:
            yield first_num + x

    def get_generator(self, first, before, middle, after):
        '''
        Checks whether it should be letters or numbers generator.
        In case of numbers checks whether the new numbering is possible
        '''
        if len(middle) < 26: 
            return self.letter_generator(first)
        n = NumberRenumerator(self.identifiers)
        if before and after:
            if int(after[0])-int(before[-1]) > len(middle): 
                return n.number_generator(first)
            else:
                raise RenumeratorError('Cannot generate numbers - not enough room for %i numbers between %s and %s.'%(len(middle), before[-1], after[0])) 
        elif before:
            return n.number_generator(first)
        else:
            return n.number_generator(first)
        
    def find_first_number(self, before, middle, after):
        """
        Returns starting number and KR what when the first number has already a letter?
        Should it generate an exception?
        """
        if not before and not after:
            return '1'
        elif before:
            return before[-1]
        elif after and before and len(middle) >= 26: 
            return str(int(after[0]) - int(before[-1])-1)
        elif after and len(middle) >= 26: 
            first = int(after[0])
            if first > len(middle):
                return str(int(after[0])-len(middle)-1)
            else:
                raise RenumeratorError('Cannot generate numbers - not enough room for %i numbers before %s.'%(len(middle), after[0])) 
        else: 
            return str(int(after[0])-1)


def renumber_section(struc, start, end, new_start):
    """
    Renumbers a part of a structure 
    from start..end, new numbers starting from new_start.
    """
    length = end - start + 1
    for old, new in zip(range(start, start + length), range(new_start, new_start+length)):
        struc.renumber_residue(str(old), str(new))
