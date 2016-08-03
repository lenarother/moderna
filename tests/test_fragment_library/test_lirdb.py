#!/usr/bin/env python
"""
Unit tests for linkers library
"""

from unittest import TestCase, main
from moderna.fragment_library.LIRdb import MakeLirFile
from moderna.fragment_library.LIR import LirRecord
from moderna.fragment_library.SearchLIR import FragmentCandidates, LirQuery, FragmentFinder
from moderna.ModernaStructure import ModernaStructure
from moderna.sequence.ModernaSequence import Sequence
from test_data import *
import os

LIR_TEST_FILE = TEST_DATA_PATH + 'lir_test_db.out'
LIR_STRUCTURES_DIR= TEST_DATA_PATH + 'lir_test_files/'
STRUC = 'pr0047Hc_with_gap_and_N.pdb'
CHAIN_LIST = TEST_DATA_PATH + 'lir_chains.txt'
# has a gap between resi 17 and 19
# has a N residue at resi 24


class LirDbTests(TestCase):

    def setUp(self):
        self.struc = ModernaStructure('file',LIR_STRUCTURES_DIR+STRUC, 'C')
        if os.path.exists(LIR_TEST_FILE):
            os.system('rm -r %s'%LIR_TEST_FILE)
            
    def tearDown(self):
        if os.path.exists(LIR_TEST_FILE):
            os.system('rm -r %s'%LIR_TEST_FILE)

    def test_get_records_from_one_chain(self):
        """Should return LIR_Records for one chain"""
        lirdb = MakeLirFile(LIR_STRUCTURES_DIR, CHAIN_LIST)
        lirdb.info[STRUC] = {}
        records = lirdb.get_records_from_one_chain(STRUC, 'C')
        self.assertEqual(len(records), 224)
        self.assertTrue(isinstance(records[0], LirRecord))
        self.assertNotEqual(str(records[0]), str(records[1]))
        
    def dont_test_nonexisting_chain(self):
        """What should happen to a chain that does not exist?"""
        records = lirdb.get_records_from_one_chain(STRUC, 'B')
        self.assertEqual(len(records), 0)
        
    def test_get_LIR_record(self):
        """Should make one LIR_Record."""
        lirdb = MakeLirFile(LIR_STRUCTURES_DIR, CHAIN_LIST)
        lirdb.info[STRUC] = {'C': {'record_number': 0,  'record_errors': 0}}
        record = lirdb.get_lir_record(self.struc['2':'5'], STRUC, 'C')
        self.assertTrue(isinstance(record, LirRecord))
        self.assertEqual(record.preceding_residue, '2')
        self.assertEqual(record.following_residue, '5')

    def test_gaps_and_unknown(self):
        """Should discard gaps and unknown residues."""
        lirdb = MakeLirFile(LIR_STRUCTURES_DIR, CHAIN_LIST) 
        lirdb.info[STRUC] = {}
        records = lirdb.get_records_from_one_chain(STRUC, 'C')
        for r in records:
            # no loops with the gap in between
            self.assertFalse(int(r.preceding_residue)<=17 and int(r.following_residue)>=19) 
            # no loops with the unknown residue in between
            self.assertFalse(int(r.preceding_residue)<=24 and int(r.following_residue)>=24)

    def test_generate_db(self):
        """Should generate a LIR database file."""
        lirdb = MakeLirFile(LIR_STRUCTURES_DIR, CHAIN_LIST)
        lirdb.generate_lir_db()
        lirdb.write_lir_db_to_file(LIR_TEST_FILE)
        self.assertTrue(os.access(LIR_TEST_FILE,os.F_OK))
        
    def test_generate_db_contents(self):
        """The generated file should be parseable."""
        # get query
        lf = FragmentFinder(self.struc['2'], self.struc['6'], Sequence('AGCU'), self.struc)
        query = lf.get_query()
        # prepare loop database
        lirdb = MakeLirFile(LIR_STRUCTURES_DIR, CHAIN_LIST)
        lirdb.generate_lir_db()
        lirdb.write_lir_db_to_file(LIR_TEST_FILE)
        # count entries in it
        lc = FragmentCandidates(query, LIR_TEST_FILE)
        lc.parse_lir_database()
        self.assertTrue(len(lc.lir_cache[1])>50)
        
    def test_get_records_from_one_structure(self):
        """Should generate a list of LIR_Record objects"""
        lirdb = MakeLirFile(LIR_STRUCTURES_DIR, CHAIN_LIST)
        records = lirdb.get_records_from_one_structure(STRUC)
        self.assertEqual(len(records),224)
        self.assertTrue(isinstance(records[0], LirRecord))
        self.assertNotEqual(str(records[0]), str(records[1]))

    
if __name__ == '__main__':
    main()
    
