#!/usr/bin/env python
"""
unit tests for different functions working with isosteric base pairs
"""


from unittest import TestCase, main
from moderna.Constants import DATA_PATH
from moderna.isosteric.Isostericity import Isostericity
from moderna.ModernaStructure import ModernaStructure

class _____IsostericityTests(TestCase):
    '''
    SWITCHED OFF 2014/08/15
    '''

    def xtest_superimpose(self):
        struct = ModernaStructure('file', 'test_data/rna_structures/1EHZ.pdb')
        iso = Isostericity((struct['14'], struct['13']), 'AC')
        rmsd = iso.rmsd       
        self.assertTrue(rmsd < iso.offset)

    def xtest_superimpose2(self):
        struct = ModernaStructure('file', 'test_data/rna_structures/1EHZ.pdb')
        iso = Isostericity((struct['20'], struct['22']), 'CU')
        rmsd = iso.rmsd       
        self.assertTrue(rmsd > iso.offset)
    
    def xtest_result_bp_numbers(self):
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
