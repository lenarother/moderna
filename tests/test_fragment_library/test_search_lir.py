#!/usr/bin/env python
#
# test_poweruser.py
#
# unit tests for base exchange functionality
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
from moderna.fragment_library.LIR import LirRecord
from moderna.ModernaStructure import ModernaStructure
from moderna.ModernaFragment import ModernaFragment
from moderna.sequence.ModernaSequence import Sequence
from moderna.fragment_library.SearchLIR import FragmentFinder, FragmentCandidates, LirQuery, LirHit, LirScoringOption
from moderna import load_model, fix_backbone, find_fragment
from test_data import *
import os
from math import pi,radians

FRAGMENT_DB = TEST_DATA_PATH + 'other/LIR_test_db'
TEST_LIR_PATH = TEST_DATA_PATH + 'lir_test_files/'

class FragmentFinderTests(TestCase):
    def setUp(self):
        self.struc = ModernaStructure('file', MINI_TEMPLATE)
        
    def tearDown(self):
        self.struc = None

    def test_init(self):
        """FragmentFinder should be initialized with a place to fit in loops"""
        lf = FragmentFinder(self.struc['2'], self.struc['6'], Sequence('AGCU'), self.struc)
        lf = FragmentFinder(self.struc['2'], self.struc['6'], Sequence('AGCU'), self.struc, 10)
        self.assertTrue(1)
        
    def test_get_query(self):
        """Should generate a LirQuery instance"""
        lf = FragmentFinder(self.struc['2'], self.struc['6'], Sequence('AGCU'), self.struc)
        q = lf.get_query()
        self.assertTrue(isinstance(q, LirQuery))
        
    def test_find_fragment_candidates(self):
        """Should return a FragmentCandidates instance"""
        lf = FragmentFinder(self.struc['2'], self.struc['6'], Sequence('AGCU'), self.struc, 10)
        cand = lf.find_fragment_candidates()
        self.assertTrue(isinstance(cand, FragmentCandidates))
        self.assertTrue(1<=len(cand)<=10)
        
    def test_find_fragment(self):
        """Should return a ModernaFragment instance with the right sequence"""
        lf = FragmentFinder(self.struc['2'], self.struc['6'], Sequence('AGCU'), self.struc, 10)
        frag = lf.find_fragment()
        frag.fix_backbone()
        self.assertTrue(isinstance(frag, ModernaFragment))
        self.assertEqual(frag.struc.get_sequence().seq_with_modifications[1:-1], 'AGCU')

    def test_find_zero_length_fragment(self):
        lf = FragmentFinder(self.struc['2'], self.struc['4'], Sequence(''), self.struc, 10)
        frag = lf.find_fragment()
        self.assertEqual(frag.struc.get_sequence(), Sequence('C_G'))
    
    def test_find_secstruc(self):
        """Secstruc can be used as an optional search criterion"""
        hairpin = ModernaStructure('file',RNA_HAIRPIN, 'D')
        lf = FragmentFinder(hairpin['30'], hairpin['40'], Sequence('CGGGCG'), self.struc, 10, secstruc="(....)")
        frag = lf.find_fragment()
        self.assertEqual(frag.struc.get_secstruc(), "((....))")

    def test_find_secstruc20(self):
        """Secstruc can be used as an optional search criterion"""
        hairpin = load_model(RNA_HAIRPIN, 'D')
        frag = find_fragment(hairpin, 29, 40, 'AUUCGUEAU', number_of_candidates=20, secstruc='(.......)')        
        self.assertEqual(frag[0].fragment_instance.struc.get_secstruc(), ".(.......).")


class FragmentCandidatesTests(TestCase):

    def setUp(self):
        self.struc = ModernaStructure('file',MINI_TEMPLATE)
        lf = FragmentFinder(self.struc['2'], self.struc['6'], Sequence('AGCU'), self.struc, lir_path=TEST_LIR_PATH)
        self.query = lf.get_query()
        self.lc = FragmentCandidates(self.query, lir_db_filename=FRAGMENT_DB)

    def tearDown(self):
        """dispose references to objects."""
        self.lc = None
        self.struc = None
        self.query = None

    def test_parse_lir_database(self):
        self.lc.parse_lir_database()
        self.assertEqual(len(self.lc.lir_cache[1]),1916)
        
    def test_create_initial_fragment_set(self):
        """Should find some hits"""
        self.assertEqual(len(self.lc), 0)
        self.lc.create_initial_fragment_set()
        self.assertTrue(len(self.lc) > 0)
        self.assertTrue(isinstance(self.lc[0], LirHit))
    
    def test_make_fast_scoring(self):
        """Should result in hits with scores"""
        self.lc.create_initial_fragment_set()
        self.assertEqual(self.lc[0].score, 0.0)
        so = LirScoringOption('fast')
        self.lc.make_fast_scoring(so, 3)
        self.assertTrue(0<=len(self.lc)<=3)
        self.assertNotEqual(self.lc[0].score, 0.0)
    
    def test_make_advanced_scoring(self):
        """Should result in hits with RMSD"""
        self.lc.create_initial_fragment_set()
        so = LirScoringOption('fast')
        self.lc.make_fast_scoring(so, 10)
        so = LirScoringOption('advanced')
        self.assertEqual(self.lc[0].rmsd, 0.0)
        self.lc.make_advanced_scoring(so)
        self.assertNotEqual(self.lc[0].rmsd, 0.0)
        
class WriteFragmentCandidatesTests(TestCase):
    #KR: load model and candidates only once --> saves lots of time
    m = load_model(TEST_DATA_PATH + 'gaps/gap7_template.pdb','D')
    candidates = m.find_fragment_candidates(m['515'], m['519'],Sequence('DDG'), 10)#, lir_path='test_data/')     
    #TODO: solve path problem: RNAModel needs to know LIR db file as optional parameter
    
    def setUp(self):
        self.remove_dir()
        
    def remove_dir(self,dirname=TEST_DATA_PATH + 'loop_candidates'):
        if os.access(TEST_DATA_PATH + 'loop_candidates', os.F_OK):
            for fn in os.listdir(dirname):
                os.remove(dirname+os.sep+fn)
            os.rmdir(dirname)
        
    def test_candidate_number(self):
        self.assertEqual(len(self.candidates), 10)
        
    def test_write_candidates_options_without_stems(self):
        """ Options for writing candidates should work 1. without stem residues"""
        self.candidates.write_fragment_candidates(TEST_DATA_PATH + 'loop_candidates', False, False, False, True)
        files = [fn for fn in os.listdir(TEST_DATA_PATH + 'loop_candidates')]
        self.assertEqual(len(files), 11)
        # check a sample for the right number of residues (3)
        m = ModernaStructure('file', TEST_DATA_PATH + 'loop_candidates/'+files[1])
        self.assertEqual(len(m), 3)
        
    def test_write_candidates_options_with_stems(self):
        """ Options for writing candidates should work 2. with stem residues"""
        self.candidates.write_fragment_candidates(TEST_DATA_PATH + 'loop_candidates', True, False, False, True)
        files = [fn for fn in os.listdir(TEST_DATA_PATH + 'loop_candidates')]
        self.assertEqual(len(files), 11)
        # check a sample for the right number of residues (5)
        m = ModernaStructure('file', TEST_DATA_PATH + 'loop_candidates/'+files[1])
        self.assertEqual(len(m), 5)
        
    def test_write_candidates_options_with_model(self):
        """ Options for writing candidates should work 3. without model"""
        self.candidates.write_fragment_candidates(TEST_DATA_PATH + 'loop_candidates', False, True, False, True)
        files = [fn for fn in os.listdir(TEST_DATA_PATH + 'loop_candidates')]
        self.assertEqual(len(files), 11)
        m = ModernaStructure('file', TEST_DATA_PATH + 'loop_candidates/'+files[1])
        self.assertEqual(len(m), 10)
        
    def _test_write_candidates_with_right_seq(self):
        """All written candidates should have the right sequence."""
        #TODO: reactivate tet for BBBuilder improvement
        self.candidates.write_fragment_candidates(TEST_DATA_PATH + 'loop_candidates', True, True, False,  False)
        i = 0 
        for fn in os.listdir(TEST_DATA_PATH + 'loop_candidates'):
            if fn.endswith('.log'): continue
            i += 1
            #if i>3: continue
            m = load_model(TEST_DATA_PATH + 'loop_candidates/'+fn)
            self.assertEqual(m.get_sequence(), Sequence('AGDDGGUUAG') )


class LirHitTests(TestCase):
    
    def setUp(self):
        self.struc = ModernaStructure('file',MINI_TEMPLATE)
        lf = FragmentFinder(self.struc['2'], self.struc['6'], Sequence('AGCU'), self.struc, lir_path=TEST_LIR_PATH)
        self.query = lf.get_query()
        self.lc = FragmentCandidates(self.query, lir_db_filename=FRAGMENT_DB)
        self.lc.create_initial_fragment_set()
        self.hit = self.lc[0]
        
    def tearDown(self):
        """dispose references to objects."""
        self.lc = None
        self.struc = None
        self.hit = None
        self.query = None
        
    def test_get_fragment(self):
        """Should return a ModernaFragment instance"""
        frag = self.hit.get_fragment()
        self.assertTrue(isinstance(frag, ModernaFragment))
        
class DihedralPerformance(TestCase):

    def setUp(self):
        self.struc = ModernaStructure('file',MINI_TEMPLATE)
        lf = FragmentFinder(self.struc['2'], self.struc['6'], Sequence('AGCU'), self.struc, lir_path=TEST_LIR_PATH)
        self.query = lf.get_query()
        self.lc = FragmentCandidates(self.query, lir_db_filename=FRAGMENT_DB)
        self.lc.create_initial_fragment_set()
        self.hit = self.lc[0]    

    def tearDown(self):
        """dispose references to objects."""
        self.lc = None
        self.struc = None
        self.hit = None
        self.query = None

    def __test_dihedral_performance(self):
        """Should calculate diff of dihedrals in radians."""
        examples = [
            (1.0,1.0,0.0),
            (0.0,0.5*pi,0.5*pi),
            (0.0,pi,pi),
            (0.0,0.5*pi*3,0.5*pi),
            (0.0,pi*2,0.0),
            (1.0,pi,pi-1),
            (1.0,0.5*pi+1,0.5*pi),
            (1.0,pi,pi-1.0),
            (1.0,0.5*pi*3,0.5*pi+1.0),
            (1.0,pi*2,1.0),
            (0.5*pi+1  ,pi,0.5*pi-1),
            (0.5*pi+1  ,0.5*pi*3,pi-1),
            (0.5*pi+1  ,pi*2,0.5*pi+1),
            (pi+1,0.5*pi*3,0.5*pi-1),
            (pi+1,pi*2,pi-1),
            (0.5*pi*3  ,pi*2,0.5*pi),
            (0.5*pi+1  ,0.5*pi+2,1),
            (radians(10), radians(220.0), radians(150)),
        ]
        for i in range(10000):
            for d1,d2,result in examples:
                self.hit.calculate_dihedrals_difference(d1,d2)
                self.hit.calculate_dihedrals_difference(d2,d1)


if __name__ == '__main__':
    main()
  
