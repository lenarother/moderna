#!/usr/bin/env python
"""
Unit tests for parsing entire linkers colection
"""

from unittest import TestCase, main
from moderna.ModernaFragment import ModernaFragment53
from moderna.fragment_library.StructureLibrary import StructureLibrary
from moderna.sequence.ModernaSequence import Sequence
from moderna.Constants import LIR_DATABASE_PATH, PATH_TO_LIR_STRUCTURES
from moderna.tests.test_data import *


class LIRDatabaseTests(TestCase):
    
    # DISABLED: slow
    def _test_all_records(self):
        """All records in the entire LIR db should be readable."""
        library = StructureLibrary(PATH_TO_LIR_STRUCTURES)
        f = open(LIR_DATABASE_PATH)
        for i,line in enumerate(f):
            print(i,line)
            if line.startswith('loop length'): 
                continue
            line=line.strip().split('\t')
            if len(line) == 17:
                st = library.get_structure_part(line[1], line[2], line[3], line[4])
                seq = Sequence(line[6])
                fr = ModernaFragment53('residues', st[line[3] : line[4]], line[2], st[line[3]], st[line[4]], None, None, None )     
                assert len(fr)-2==int(line[0])
                fr_seq = fr.get_sequence()
                print('1',seq,'2',fr_seq)
                assert seq==fr_seq
            if i % 100==0: 
                print(i)


if __name__ == '__main__':
    main()


