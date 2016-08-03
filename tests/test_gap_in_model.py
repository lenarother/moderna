#!/usr/bin/env python
"""
Unit Tests for linker insertion and backbone continuity
"""

from unittest import main, TestCase
from moderna.ModernaStructure import ModernaStructure
from moderna.fragment_library.SearchLIR import *
from moderna import *
from moderna.ModernaSequence import Sequence
import os

FRAGMENTS_NUMBER = 5000


class GapInModelTests(TestCase):
    def test_gap1_model(self):
        """
        Should create model with no gaps.
        Template: 1B23 (fragment)
        Target: 1EHZ (fragment)
        """
        t = load_template('test_data/gaps/gap1_template.pdb','B')
        a = load_alignment('test_data/gaps/gap1_alignment.fasta')
        m = create_model(t, a)
        os.system('mkdir -p test_data/gaps/gap1_candidates')
        old_candidates = os.listdir('test_data/gaps/gap1_candidates')
        for x in old_candidates: 
            if x.startswith('model'): os.remove('test_data/gaps/gap1_candidates/'+x)
        m.write_pdb_file('test_data/gaps/gap1_candidates/model1.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('AGDDGGG'))
        

#    def test_gap1_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 1B23 (fragment)
#        Target: 1EHZ (fragment)
#        """
#        m = load_model('test_data/gaps/gap1_template.pdb','B')
#        candidates = m.find_loop_candidates(m['915'], m['919'],Sequence('DDG'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap1_candidates')
#        old_candidates = os.listdir('test_data/gaps/gap1_candidates')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap1_candidates/'+x)
#        candidates.write_loop_candidates('test_data/gaps/gap1_candidates',True, True,  True)
#        self.assertEqual(len(candidates), FRAGMENTS_NUMBER)


    def test_gap2_model(self):
        """
        Should create model with no gaps.
        Template: 1H4S  T (fragment)
        Target: 1EHZ A (fragment)
        """
        t = load_template('test_data/gaps/gap2_template.pdb','T')
        a = load_alignment('test_data/gaps/gap2_alignment.fasta')
        m = create_model(t, a)
        os.system('mkdir -p test_data/gaps/gap2_candidates')
        old_candidates = os.listdir('test_data/gaps/gap2_candidates')
        for x in old_candidates: 
            if x.startswith('model'): os.remove('test_data/gaps/gap2_candidates/'+x)
        write_model(m,'test_data/gaps/gap2_candidates/model2.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('CAGAAUUC'))
    
    
#    def test_gap2_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 1H4S  T (fragment)
#        Target: 1EHZ A (fragment)
#        """
#        m = load_model('test_data/gaps/gap2_template.pdb','T')
#        candidates = m.find_loop_candidates(m['64'], m['67'],Sequence('GAA'),FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap2_candidates')
#        old_candidates = os.listdir('test_data/gaps/gap2_candidates')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap2_candidates/'+x)
#        candidates.write_loop_candidates('test_data/gaps/gap2_candidates', True, True,  True)
#        self.assertEqual(len(candidates), FRAGMENTS_NUMBER)

  
    def test_gap3_model(self):
        """
        Should create model with no gaps.
        Template: 1QF6 B  (fragment)
        Target: 1B23 R  (fragment)
        """
        t = load_template('test_data/gaps/gap3_template.pdb','B')
        a = load_alignment('test_data/gaps/gap3_alignment.fasta')
        m = create_model(t, a)
        os.system('mkdir -p test_data/gaps/gap3_candidates')
        old_candidates = os.listdir('test_data/gaps/gap3_candidates')
        for x in old_candidates:
           if x.startswith('model'): os.remove('test_data/gaps/gap3_candidates/'+x)
        write_model(m,'test_data/gaps/gap3_candidates/model3.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('CUAGUC'))
    
    
#    def test_gap3_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 1QF6 B  (fragment)
#        Target: 1B23 R  (fragment)
#        """
#        m = load_model('test_data/gaps/gap3_template.pdb','B')
#        candidates = m.find_loop_candidates(m['45'], m['49'],Sequence('AG'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap3_candidates')
#        old_candidates = os.listdir('test_data/gaps/gap3_candidates')
#        for x in old_candidates: 
#            if x.startwith('candidate'): os.remove('test_data/gaps/gap3_candidates/'+x)
#        candidates.write_loop_candidates('test_data/gaps/gap3_candidates', True, True,  True)
#        self.assertEqual(len(candidates), FRAGMENTS_NUMBER)
    
    
    def test_gap4_model(self):
        """
        Should create model with no gaps.
        Template: 1WZ2 C  (fragment)
        Target: 1B23 R  (fragment)
        """
        t = load_template('test_data/gaps/gap4_template.pdb','C')
        a = load_alignment('test_data/gaps/gap4_alignment.fasta')
        m = create_model(t, a)
        os.system('mkdir -p test_data/gaps/gap4_candidates1')
        old_candidates = os.listdir('test_data/gaps/gap4_candidates1')
        for x in old_candidates: 
            if x.startswith('model'): os.remove('test_data/gaps/gap4_candidates1/'+x)
        write_model(m,'test_data/gaps/gap4_candidates1/model4.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('AGCGGDDAU'))  


#    def test_gap4_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 1WZ2 C  (fragment)
#        Target: 1B23 R  (fragment)
#        """
#        m = load_model('test_data/gaps/gap4_template.pdb','C')
#        candidates1 = m.find_loop_candidates(m['915'], m['920'],Sequence('CG'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap4_candidates1')
#        old_candidates = os.listdir('test_data/gaps/gap4_candidates1')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap4_candidates1/'+x)
#        candidates1.write_loop_candidates('test_data/gaps/gap4_candidates1', True, True,  True)
#        candidates2 = m.find_loop_candidates(m['920'], m['925'],Sequence('DD'),FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap4_candidates2')
#        old_candidates = os.listdir('test_data/gaps/gap4_candidates2')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap4_candidates2/'+x)
#        candidates2.write_loop_candidates('test_data/gaps/gap4_candidates2', True, True,  True)
#        self.assertEqual(len(candidates1), FRAGMENTS_NUMBER)
#        self.assertEqual(len(candidates2), FRAGMENTS_NUMBER)
    

    def test_gap5_model(self):
        """
        Should create model with no gaps.
        Template: 1WZ2 C  (fragment)
        Target: 1B23 R  (fragment)
        """
        # Won't find gap automaticly.
        t = load_template('test_data/gaps/gap5_template.pdb','C')
        a = load_alignment('test_data/gaps/gap5_alignment.fasta')
        m = create_model(t, a)
        os.system('mkdir -p test_data/gaps/gap5_candidates')
        old_candidates = os.listdir('test_data/gaps/gap5_candidates')
        for x in old_candidates: 
            if x.startswith('model'): os.remove('test_data/gaps/gap5_candidates/'+x)
        write_model(m,'test_data/gaps/gap5_candidates/model5.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('GUCUAGU'))
    
    
#    def test_gap5_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 1WZ2 C  (fragment)
#        Target: 1B23 R  (fragment)
#        """
#        # Won't find gap automaticly.
#        m = load_model('test_data/gaps/gap5_template.pdb','C')
#        candidates = m.find_loop_candidates(m['946'], m['960'],Sequence('CUA'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap5_candidates')
#        old_candidates = os.listdir('test_data/gaps/gap5_candidates')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap5_candidates/'+x)
#        candidates.write_loop_candidates('test_data/gaps/gap5_candidates',True, True,  True)
#        self.assertEqual(len(candidates), FRAGMENTS_NUMBER)

    
    
    def test_gap6_model(self):
        """
        Should create model with no gaps.
        Template: 2V0G B  (fragment)
        Target: 1GTS B  (fragment)
        """
        # Template has letters in original numeration.
        t = load_template('test_data/gaps/gap6_template.pdb','B')
        a = load_alignment('test_data/gaps/gap6_alignment.fasta')
        m = create_model(t, a)    
        os.system('mkdir -p test_data/gaps/gap6_candidates')
        old_candidates = os.listdir('test_data/gaps/gap6_candidates')
        for x in old_candidates: 
            if x.startswith('model'): os.remove('test_data/gaps/gap6_candidates/'+x)
        write_model(m,'test_data/gaps/gap6_candidates/model6.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('UUCCGA'))   
    
    
#    def test_gap6_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 2V0G B  (fragment)
#        Target: 1GTS B  (fragment)
#        """
#        # Template has letters in original numeration.
#        m = load_model('test_data/gaps/gap6_template.pdb','B')
#        #candidates = m.find_loop_candidates(m['47E'], m['50'],Sequence('CC'), FRAGMENTS_NUMBER)
#        candidates = m.find_loop_candidates(m['2'], m['6'],Sequence('CC'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap6_candidates')
#        old_candidates = os.listdir('test_data/gaps/gap6_candidates')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap6_candidates/'+x)
#        candidates.write_loop_candidates('test_data/gaps/gap6_candidates', True, True,  True)
#        self.assertEqual(len(candidates), FRAGMENTS_NUMBER)
    

    def test_gap7_model(self):
        """
        Should create model with no gaps.
        Template: 2CV2 D  (fragment)
        Target: 1WZ2 C  (fragment)
        """
        t = load_template('test_data/gaps/gap7_template.pdb','D')
        a = load_alignment('test_data/gaps/gap7_alignment.fasta')
        m = create_model(t, a)
        os.system('mkdir -p test_data/gaps/gap7_candidates1')
        old_candidates = os.listdir('test_data/gaps/gap7_candidates1')
        for x in old_candidates:
           if x.startswith('model'): os.remove('test_data/gaps/gap7_candidates1/'+x)
        write_model(m,'test_data/gaps/gap7_candidates1/model7.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('AGCCUGGUCAAA')) 
        
#    def test_gap7_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 2CV2 D  (fragment)
#        Target: 1WZ2 C  (fragment)
#        """
#        m = load_model('test_data/gaps/gap7_template.pdb','D')
#        #candidates1 = m.find_loop_candidates(m['515'], m['519'],Sequence('CCUG'), FRAGMENTS_NUMBER)
#        candidates1 = m.find_loop_candidates(m['2'], m['5'],Sequence('CCUG'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap7_candidates1')
#        old_candidates = os.listdir('test_data/gaps/gap7_candidates1')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap7_candidates1/'+x)
#        candidates1.write_loop_candidates('test_data/gaps/gap7_candidates1', True, True,  True)
#        #candidates2 = m.find_loop_candidates(m['519'], m['521'],Sequence('UCA'), FRAGMENTS_NUMBER)
#        candidates2 = m.find_loop_candidates(m['5'], m['8'],Sequence('UCA'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap7_candidates2')
#        old_candidates = os.listdir('test_data/gaps/gap7_candidates2')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap7_candidates2/'+x)
#        candidates2.write_loop_candidates('test_data/gaps/gap7_candidates2', True, True,  True)
#        #candidates3 = m.find_loop_candidates(m['515'], m['521'],Sequence('CCUGGUCA'), FRAGMENTS_NUMBER)
#        candidates3 = m.find_loop_candidates(m['2'], m['8'],Sequence('CCUGGUCA'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap7_candidates3')
#        old_candidates = os.listdir('test_data/gaps/gap7_candidates3')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap7_candidates3/'+x)
#        candidates3.write_loop_candidates('test_data/gaps/gap7_candidates3', True, True,  True)
#        self.assertEqual(len(candidates1), FRAGMENTS_NUMBER)
#        self.assertEqual(len(candidates2), FRAGMENTS_NUMBER)
#        self.assertEqual(len(candidates3), FRAGMENTS_NUMBER)
    
    
    def test_gap8_model(self):
        """
        Should create model with no gaps.
        Template: 2BYT B  (fragment)
        Target: 2QNH z (fragment)
        """
        t = load_template('test_data/gaps/gap8_template.pdb','B')
        a = load_alignment('test_data/gaps/gap8_alignment.fasta')
        m = create_model(t, a)
        os.system('mkdir -p test_data/gaps/gap8_candidates')
        old_candidates = os.listdir('test_data/gaps/gap8_candidates')
        for x in old_candidates: 
            if x.startswith('model'): os.remove('test_data/gaps/gap8_candidates/'+x)
        write_model(m,'test_data/gaps/gap8_candidates/model8.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('GCAGCCU'))
        
        
#    def test_gap8_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 2BYT B  (fragment)
#        Target: 2QNH z (fragment)
#        """
#        m = load_model('test_data/gaps/gap8_template.pdb','B')
#        candidates = m.find_loop_candidates(m['13'], m['16'],Sequence('AGC'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/gap8_candidates1')
#        old_candidates = os.listdir('test_data/gaps/gap8_candidates')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/gap8_candidates/'+x)
#        candidates.write_loop_candidates('test_data/gaps/gap8_candidates', True, True,  True)
#        self.assertEqual(len(candidates), FRAGMENTS_NUMBER)
    
    
    def test_gap_structure1_model(self):
        """
        Should create model with no gaps.
        Template: 1H4S T  (whole structure)
        Target: 1EHZ A (whole structure)
        """
        t = load_template('test_data/gaps/1h4s_T.pdb','T')
        a = load_alignment('test_data/gaps/structure1_alignment.fasta')
        m = create_model(t, a)
        os.system('mkdir -p test_data/gaps/structure1_candidates1')
        old_candidates = os.listdir('test_data/gaps/structure1_candidates1')
        for x in old_candidates:
            if x.startswith('model'): os.remove('test_data/gaps/structure1_candidates1/'+x)
        write_model(m,'test_data/gaps/structure1_candidates1/model_structure1.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('GCGGAUUUALCUCAGDDGGGAGAGCRCCAGABU#AAYAP?UGGAG7UC?UGUGTPCG"UCCACAGAAUUCGCACCA'))
    
#    def test_gap_structure1_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 1H4S T  (whole structure)
#        Target: 1EHZ A (whole structure)
#        """
#        m = load_model('test_data/gaps/1h4s_T.pdb','T')
#        #print [r.identifier for r in m]           
#        #candidates1 = m.find_loop_candidates(m['17'], m['19'],Sequence('DDG'),FRAGMENTS_NUMBER)
#        candidates1 = m.find_loop_candidates(m['14'], m['17'],Sequence('DDG'),FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/structure1_candidates1')
#        old_candidates = os.listdir('test_data/gaps/structure1_candidates1')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/structure1_candidates1/'+x)
#        candidates1.write_loop_candidates('test_data/gaps/structure1_candidates1', True, True,  True)
#        #candidates2 = m.find_loop_candidates(m['19'], m['23'],Sequence('GA'),FRAGMENTS_NUMBER)
#        candidates2 = m.find_loop_candidates(m['17'], m['21'],Sequence('GA'),FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/structure1_candidates2')
#        old_candidates = os.listdir('test_data/gaps/structure1_candidates2')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/structure1_candidates2/'+x)
#        candidates2.write_loop_candidates('test_data/gaps/structure1_candidates2', True, True,  True)
#        #candidates3 = m.find_loop_candidates(m['46'], m['49'],Sequence('7U'), FRAGMENTS_NUMBER)
#        candidates3 = m.find_loop_candidates(m['42'], m['45A'],Sequence('7U'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/structure1_candidates3')
#        old_candidates = os.listdir('test_data/gaps/structure1_candidates3')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/structure1_candidates3/'+x)
#        candidates3.write_loop_candidates('test_data/gaps/structure1_candidates3', True, True,  True)
#        #candidates4 = m.find_loop_candidates(m['64'], m['67'],Sequence('GAA'), FRAGMENTS_NUMBER)
#        candidates4 = m.find_loop_candidates(m['60'], m['63A'],Sequence('GAA'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/structure1_candidates4')
#        old_candidates = os.listdir('test_data/gaps/structure1_candidates4')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/structure1_candidates4/'+x)
#        candidates4.write_loop_candidates('test_data/gaps/structure1_candidates4', True, True,  True)
#        self.assertEqual(len(candidates1), FRAGMENTS_NUMBER)
#        self.assertEqual(len(candidates2), FRAGMENTS_NUMBER)
#        self.assertEqual(len(candidates3), FRAGMENTS_NUMBER)
#        self.assertEqual(len(candidates4), FRAGMENTS_NUMBER)
    
    
    def test_gap_structure2_model(self):
        """
        Should create model with no gaps.
        Template: 1H3E B  (whole structure)
        Target: 1WZ2 C (whole structure)
        """
        t = load_template('test_data/gaps/1h3e_B.pdb','B')
        print(examine_structure(t))
        print(t['16'])
        a = load_alignment('test_data/gaps/structure2_alignment.fasta')
        print(a)
        m = create_model(t, a)
        os.system('mkdir -p test_data/gaps/structure2_candidates1')
        old_candidates = os.listdir('test_data/gaps/structure2_candidates1')
        for x in old_candidates: 
            if x.startswith('model'): os.remove('test_data/gaps/structure2_candidates1/'+x)
        print(m)
        write_model(m,'test_data/gaps/structure2_candidates1/model_structure2.pdb')
        self.assertEqual(m.get_sequence(),  Sequence('GCGGGGGUUGCCGAGCCUGGUCAAAGGCGGGGGACUCAAGAUCCCCUCCCGUA.GGGUUCCGGGGUUCGAAUCCCCGCCCCCGCACCA'))
    
    
#    def test_gap_structure2_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 1H3E B  (whole structure)
#        Target: 1WZ2 C (whole structure)
#        """
#        m = load_model('test_data/gaps/1h3e_B.pdb','B')
#        candidates1 = m.find_loop_candidates(m['15'], m['19'],Sequence('CCUG'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/structure2_candidates1')
#        old_candidates = os.listdir('test_data/gaps/structure2_candidates1')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/structure2_candidates1/'+x)
#        candidates1.write_loop_candidates('test_data/gaps/structure2_candidates1', True, True,  True)
#        candidates2 = m.find_loop_candidates(m['65'], m['68'],Sequence('CCC'), FRAGMENTS_NUMBER)
#        os.system('mkdir -p test_data/gaps/structure2_candidates2')
#        old_candidates = os.listdir('test_data/gaps/structure2_candidates2')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/structure2_candidates2/'+x)
#        candidates2.write_loop_candidates('test_data/gaps/structure2_candidates2', True, True,  True)
#        self.assertEqual(len(candidates1), FRAGMENTS_NUMBER)
#        self.assertEqual(len(candidates2), FRAGMENTS_NUMBER)
    

#    def test_gap_structure3_model(self):
#        """
#        Should create model with no gaps.
#        Template: 1QF6 B  (whole structure)
#        Target:  (whole structure)
#        """
#        t = load_template('/home/lenam/repos/model_verification/examples/trna2008/templates/1qf6_B.pdb','B')
#        print t.get_sequence()
#        a = load_alignment('test_data/gaps/structure3_alignment.fasta')
#        print a
#        m = create_model(t, a)
#        print m.get_sequence()
#        print m
#        os.system('mkdir -p test_data/gaps/structure3_candidates')
#        old_candidates = os.listdir('test_data/gaps/structure3_candidates')
#        for x in old_candidates: 
#            if x.startswith('model'): os.remove('test_data/gaps/structure3_candidates/'+x)
#        write_model(m,'test_data/gaps/structure3_candidates/model_structure3.pdb')
#        self.assertEqual(m.get_sequence(),  Sequence('GGAGCGG4AGUUCAGDCGGDDAGAAUACCUGCCUQUC/CGCAGGGG7UCGCGGGTPCGAGUCCCGPCCGUUCCGCCA'))
#
#    def test_gap_structure3_candidates(self):
#        """
#        Should create model with no gaps.
#        Template: 1QF6 B  (whole structure)
#        Target:  (whole structure)
#        """
#        m = load_model('test_data/gaps/1qf6_B.pdb','B')
#        #import moderna.ClashRecognizer
#        #cr = ClashRecognizer()
#        #print cr.find_clashes_in_residues(m)
#        candidates1 = m.find_loop_candidates(m['19'], m['22'],Sequence('DDA'), FRAGMENTS_NUMBER)
#        #for x in candidates1: print x
#        os.system('mkdir -p test_data/gaps/structure3_candidates')
#        old_candidates = os.listdir('test_data/gaps/structure3_candidates')
#        for x in old_candidates: 
#            if x.startswith('candidate'): os.remove('test_data/gaps/structure3_candidates/'+x)
#        candidates1.write_loop_candidates('test_data/gaps/structure3_candidates', True, True,  True)
#        self.assertEqual(len(candidates1), FRAGMENTS_NUMBER)
#
#
#    def test_gap_structure4_model(self):
#        """
#        Should create Lukaszs model with no gaps.
#        Should add harpin to the 3ftf helis.
#        """
#        t = load_template('test_data/gaps/3ftf_A_Lukasz_template.pdb','A')
#        a = load_alignment('test_data/gaps/structure4_Lukasz_alignment.fasta')
#        m = create_model(t, a)
#        os.system('mkdir -p test_data/gaps/structure4_candidates')
#        old_candidates = os.listdir('test_data/gaps/structure4_candidates')
#        for x in old_candidates: 
#            if x.startswith('model'): os.remove('test_data/gaps/structure4_candidates/'+x)
#        write_model(m,'test_data/gaps/structure4_candidates/model_structure4.pdb')
#        self.assertEqual(m.get_sequence(),  Sequence('CGCGACGGACGGAAAGACCCCUAUCCGUCGCG'))
        
        
if __name__ == '__main__':
    main()
  
