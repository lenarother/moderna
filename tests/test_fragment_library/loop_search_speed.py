#!/usr/bin/env python
"""
Unit Tests for linker searching speed
"""

# Search for fragment for loop in 1ehz
# Stem residues: 3,8   

from moderna.Constants import DATA_PATH
from moderna.ModernaStructure import *
from moderna.fragment_library.LIR import LIR
from moderna.fragment_library.SearchLIR import *
from moderna.RNAModel import RNAModel 
from moderna.tests.test_data import TEST_DATA_PATH
import cProfile

NUMBER_OF_LOOP_CANDIDATES = 20

def check_lir():
    """Searches for fragments fitting into 1ehz."""
    st=ModernaStructure('file',TEST_DATA_PATH + 'pdb1ehz.ent')

    l=LIR()
    l.set_residues(st[15], st[20])
    l.check_presence_of_atoms()
    l.set_vectors()

    record=LIR_Record(
        loop_length = 4,
        structure = TEST_DATA_PATH + 'pdb1ehz.ent',
        chain = 'A',
        preceding_resi = 15,
        sequence = ['G','A','U','U'],
        x=l.get_x_value(),
        y=l.get_y_value(),
        dist_stem = l.get_dist_stem(),
        beta = l.get_beta(),
        dihedral = l.get_dihedral()
        )
                        
    s=LoopSearch(record,TEST_DATA_PATH + 'LIR_test')
    # loops = s.get_loop_candidates(NUMBER_OF_LOOP_CANDIDATES)
    loops = s.get_loop_candidates(NUMBER_OF_LOOP_CANDIDATES)

    print loops
    for loop in loops: print loop


def check_insert_loop():
    m = RNAModel(model_chain_name='A', data_type='file', data=TEST_DATA_PATH + 'pdb1ehz.ent')
    m.sarch_loop_fragment(m[3], m[8], sequence=['A','A','A','A'], number_of_candidates=5)


cProfile.run("check_lir()")
#check_lir()
#check_insert_loop()
