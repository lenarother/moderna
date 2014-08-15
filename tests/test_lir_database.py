#!/usr/bin/env python
#
# test_lir_database.py
#
# unit tests for parsing entire linkers colection
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

from unittest import TestCase, main
from moderna.ModernaStructure import ModernaStructure
from moderna.ModernaFragment import ModernaFragment53
from moderna.StructureLibrary import StructureLibrary
from moderna.ModernaSequence import Sequence
from moderna.Constants import LIR_DATABASE_PATH, PATH_TO_LIR_STRUCTURES
from test_data import *

class LIRDatabaseTests(TestCase):
    
    def test_all_records(self):
        """All records in the entire LIR db should be readable."""
        library = StructureLibrary(PATH_TO_LIR_STRUCTURES)
        f = open(LIR_DATABASE_PATH)
        for i,line in enumerate(f):
            print i,line
            if line.startswith('loop length'): continue
            line=line.strip().split('\t')
            if len(line) == 17:
                st = library.get_structure_part(line[1], line[2], line[3], line[4])
                seq = Sequence(line[6])
                print 'SEQ', seq
                fr = ModernaFragment53('residues', st[line[3] : line[4]], line[2], st[line[3]], st[line[4]], None, None, None )     
                assert len(fr)-2 == int(line[0])
                fr_seq = fr.get_sequence()
                print '1', seq, '2', fr_seq
                assert seq == fr_seq
            if i%100==0: 
                print i
                

if __name__ == '__main__':
    main()


