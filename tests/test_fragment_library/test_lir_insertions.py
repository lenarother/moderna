#!/usr/bin/env python
"""
Unit tests for insertion of linkers from library

tests whether LIR finds its own fragments
"""

from unittest import TestCase, main
from moderna.ModernaStructure import ModernaStructure
from moderna.fragment_library.SearchLIR import FragmentFinder
from moderna.sequence.ModernaSequence import Sequence
from random import random
from moderna.Constants import DATA_PATH, PATH_TO_LIR_STRUCTURES, SINGLE_STRAND
from moderna import *
from moderna.ModernaSuperimposer import ModernaSuperimposer

FRAGMENTS_TO_INSERT = [
    (0 ,'pr0010Hb.pdb','B','116','117',''), 
    (5 ,'pr0010Hb.pdb','B','107','113','AUGGU') , 
    (6 ,'pr0010Hb.pdb','B','145','152','CCAUUG'), 
    (7 ,'ur0052H.pdb','X','21','29','UAAUCGC'), 
    (12,'ur0052H.pdb','X','30','43','GGAUAUGGCACG')
]


class LIRInsertionTests(TestCase):
    
    def test_insert_in_itself(self):
        """Best fragment for LIR structures should come from itself."""
        for length, struc, chain, stem5, stem3, seq in FRAGMENTS_TO_INSERT:
            m = load_model(PATH_TO_LIR_STRUCTURES+struc, chain)
            lf = FragmentFinder(m[stem5], m[stem3], Sequence(seq), m, 2)
            candidates = lf.find_fragment_candidates()
            if len(candidates)>0:
                best = candidates[0]
                self.assertEqual(best.structure, struc)
                self.assertEqual(best.chain, chain)
                self.assertEqual(best.preceding_residue, stem5)
                self.assertEqual(best.following_residue, stem3)
                self.assertAlmostEqual(best.rmsd, 0.000, 3)
                
    def get_atoms(self, mol):
        """Returns a list of atoms for superposition"""
        result = []
        atoms = ["P", "C5'", "N*", "C1'", "C2'", "C5", "C4'", "OP1", "OP2"]
        for r in mol:
            for a in atoms:
                result.append(r[a])
        return result
        
    def test_stems_same_position(self):
        """Stem residue ends should stay where they are"""
        MAX = 10
        for i in range(MAX):
            helix = load_model(SINGLE_STRAND, 'A')
            m = load_model(SINGLE_STRAND, 'A')
            lc = find_fragment(m, '1475', '1485', Sequence('GCCCGUGAC'), MAX)
            cand = lc[i]
            insert_fragment(m, cand)
            dist1=helix['1475']["C5'"] - m['1475']["C5'"]
            dist2=helix['1485']["C3'"] - m['1485']["C3'"]
            self.assertAlmostEqual(dist1, 0.0, 3)
            self.assertAlmostEqual(dist2, 0.0, 3)
        
        
    def test_insert_in_helix(self):
        """Evaluates how a helical single strand is handled by LIR"""
        helix = load_model(SINGLE_STRAND, 'A')
        m = load_model(SINGLE_STRAND, 'A')
        lc = find_fragment(m, '1475', '1485', Sequence('GCCCGUGAC'), 10)
        
        # calculate all-atom-rmsds
        result = []
        i =1
        for cand in lc:
            insert_fragment(m, cand)
            fixed = self.get_atoms(helix)
            moved = self.get_atoms(m)
            s = ModernaSuperimposer(fixed, moved, m.get_all_atoms())
            result.append(s.rmsd)
            i+=1

        self.assertAlmostEqual(result[0], 0.000, 3)
        avg_rms = sum(result)/(1.0*len(result))
        self.assertTrue(avg_rms < 3.25)
           
if __name__ == '__main__':
    main()
    
