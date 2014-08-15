#!/usr/bin/env python
#
# test_isostericity.py
#
# unit tests for different functions working with isosteric base pairs
#
__author__ = "Pawel Skiba, Magdalena Rother, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Pawel Skiba"
__email__ = "pw.skiba@gmail.com"
__status__ = "Prototype"

from unittest import TestCase, main
from moderna.Constants import DATA_PATH
from moderna.isosteric.Isostericity import Isostericity
from moderna.ModernaStructure import ModernaStructure

class _____IsostericityTests(TestCase):
    '''
    SWITCHED OFF 2014/08/15
    '''

    def test_superimpose(self):
        struct = ModernaStructure('file', 'test_data/rna_structures/1EHZ.pdb')
        iso = Isostericity((struct['14'], struct['13']), 'AC')
        rmsd = iso.rmsd       
        self.assertTrue(rmsd < iso.offset)

    def test_superimpose2(self):
        struct = ModernaStructure('file', 'test_data/rna_structures/1EHZ.pdb')
        iso = Isostericity((struct['20'], struct['22']), 'CU')
        rmsd = iso.rmsd       
        self.assertTrue(rmsd > iso.offset)
    
    def test_result_bp_numbers(self):
        struct = ModernaStructure('file', 'test_data/rna_structures/1EHZ.pdb')
        iso = Isostericity((struct['12'], struct['23']), 'CG')  
        source = [12, 23]
        result = [iso.result_bp[0].number, iso.result_bp[1].number]
        self.assertEqual(source,result)

    #TODO: Superposition of a near isosteric pair.
    #TODO: Superposition of a non-isosteric pair.
    #TODO: testing something that causes an exception.
    

if __name__ == '__main__':
    main()
