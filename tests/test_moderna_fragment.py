#!/usr/bin/env python
#
# test_moderna_fragment.py
#
# unit tests for ModernaStructure.ModernaFragment
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
from moderna.ModernaFragment import ModernaFragment,  ModernaFragment5,  \
    ModernaFragment3, ModernaFragment53, AnchorResidue, AnchorBuildRule, \
    keep_nothing, keep_first_last
from moderna.RNAModel import RnaModel
from moderna.ModernaStructure import ModernaStructure
from moderna.ModernaSequence import Sequence
from moderna.Errors import ModernaFragmentError, ModernaSuperimposerError
from moderna.Constants import HELIX
from moderna import load_template, create_model, copy_some_residues, create_fragment
from test_data import *
import re


TEST_RULE = AnchorBuildRule('model', ("C3'", "O3'"), ("O2'", "C2'"))

class AnchorResidueTests(TestCase):
    
    def setUp(self):
        self.struc = ModernaStructure('file', MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
        self.anchor = AnchorResidue(self.struc['1'], self.struc['3'], 
                      superpos_atoms=("C4'", "P", "N*"), build_rule=TEST_RULE)
                      
    def test_attributes(self):
        """Create an AnchorResidue"""
        self.assertEqual(self.anchor.fixed_resi, self.struc['1'])
        self.assertEqual(self.anchor.mobile_resi, self.struc['3'])
        self.assertEqual(len(list(self.anchor.fixed_superposition_atoms)), 3)
        self.assertEqual(len(list(self.anchor.mobile_superposition_atoms)), 3)
        for atom in self.anchor.fixed_superposition_atoms:
            self.assertEqual(atom.parent, self.anchor.fixed_resi)
        for atom in self.anchor.mobile_superposition_atoms:
            self.assertEqual(atom.parent, self.anchor.mobile_resi)
        
    def test_validate(self):
        """improper parameters raise errors."""
        self.assertRaises(ModernaFragmentError, AnchorResidue, \
                    None, self.struc['3'], ("C4'", "P", "N*"), TEST_RULE)
        self.assertRaises(ModernaFragmentError, AnchorResidue, \
                    self.struc['1'], None, ("C4'", "P", "N*"), TEST_RULE)
        self.assertRaises(ModernaFragmentError, AnchorBuildRule, \
                    'template', ("C3'", "O3'"), ("O2'", "C2'"))
        
    def test_get_anchored_residue(self):
        """Prepares a new residue object."""
        fixed_before = self.struc['1'].child_list[:]
        mobile_before = self.struc['3'].child_list[:]
        prep = self.anchor.get_anchored_residue()
        self.assertNotEqual(prep, self.struc['1'])
        self.assertNotEqual(prep, self.struc['3'])
        self.assertEqual(prep.short_abbrev, self.struc['1'].short_abbrev)
        self.assertEqual(fixed_before, self.struc['1'].child_list)
        self.assertEqual(mobile_before, self.struc['3'].child_list)
        self.assertEqual(str(prep["C2'"].coord), str(self.struc['3']["C2'"].coord))
        self.assertEqual(str(prep["C3'"].coord), str(self.struc['1']["C3'"].coord))
        
        
        
class ModernaFragmentTests(TestCase):

    def setUp(self):
        self.s1 = ModernaStructure('file',SMALL_FRAGMENT, seq=Sequence("GCGG"))
        self.s2 = ModernaStructure('file',SMALL_FRAGMENT, seq=Sequence("GCGG"))
        self.s3 = ModernaStructure('file',MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
        
    def test_init(self):
        """Should have ModernaStructure as an attribute."""
        f1 = ModernaFragment(self.s1)
        self.assertEqual(f1.struc.get_sequence(), Sequence('GCGG'))
        
    def test_attributes(self):
        """object attributes should have been set up correctly."""
        f1 = ModernaFragment(self.s1)
        f2 = ModernaFragment(self.s2, new_sequence=Sequence('GCLG'))
        modfrag = ModernaFragment(self.s3, new_sequence=Sequence('GL7L7L7ULL7L7AG'))
        self.assertEqual(len(modfrag.struc), 15)
        self.assertEqual(f1.struc.get_sequence(),Sequence("GCGG"))
        self.assertEqual(f2.struc.get_sequence(),Sequence("GCGG"))
        self.assertTrue(str(f1))
        
    def test_anchor_residues(self):
        """Abstract superclass should have no anchors yet."""
        f1 = ModernaFragment(self.s1)
        self.assertEqual(f1.anchor_residues, [])

    def test_nonanchor_residues(self):
        """Stems defined in subclasses, thus should return all residues."""
        f1 = ModernaFragment(self.s1)
        self.assertEqual(len(f1.nonanchor_residues), 4)
        
    def test_superimpose(self):
        """Should just exist as an abstract method."""
        f1 = ModernaFragment(self.s1)
        self.assertRaises(ModernaSuperimposerError, f1.superimpose)
        
    def test_prepare_anchor_residues(self):
        """Should exist as an abstract method."""
        f1 = ModernaFragment(self.s1)        
        self.assertRaises(ModernaSuperimposerError, f1.prepare_anchor_residues)
    
    def test_renumber(self):
        """Should exist as an abstract method."""
        f1 = ModernaFragment(self.s1)
        f1.renumber()
        self.assertEqual([r.identifier for r in self.s1], ['1', '2', '3', '4'])

    def test_abstract_methods(self):
        """Some placeholder methods should be there."""
        f1 = ModernaFragment(self.s1)
        f1.get_resi_to_remove(self.s2)
        f1.fix_backbone()
    
    def test_has_clashes(self):
        f1 = ModernaFragment(self.s1)
        clashes = f1.has_clashes(list(self.s2))
        self.assertEqual(len(clashes), 8)
        struc = ModernaStructure('file',FRAGMENT1)
        clashes = f1.has_clashes(list(struc))
        self.assertEqual(len(clashes), 0)
        
    def test_apply_sequence(self):
        """Should change the sequence of the entire fragment."""
        frag = ModernaFragment(self.s1, new_sequence=Sequence('AAAA'))
        frag.apply_seq()
        self.assertEqual(self.s1.get_sequence(), Sequence('AAAA'))
        # test lengths that do not fit
        frag = ModernaFragment(self.s2, new_sequence=Sequence('AAA'))
        self.assertRaises(ModernaFragmentError,frag.apply_seq)
        frag = ModernaFragment(self.s2, new_sequence=Sequence('AAAAA'))
        self.assertRaises(ModernaFragmentError,frag.apply_seq)

    def test_apply_sequence_modified(self):
        """Applying a sequence with modifications should work."""
        f2 = ModernaFragment(self.s2, new_sequence=Sequence('GCLG'))
        f2.apply_seq()
        self.assertEqual(self.s2.get_sequence(),Sequence("GCLG"))
        modfrag = ModernaFragment(self.s3, new_sequence=Sequence('GL7L7L7ULL7L7AG'))
        modfrag.apply_seq()
        self.assertEqual(self.s3.get_sequence(), Sequence('GL7L7L7ULL7L7AG'))
        
    def test_broken_fragment_fail(self):
        """Chain breaks and unknown residues cause errors."""
        broken = ModernaStructure('file', BROKEN_FRAGMENT)
        self.assertRaises(ModernaFragmentError, ModernaFragment, broken)
        unknown = ModernaStructure('file', UNKNOWN_FRAGMENT)
        self.assertRaises(ModernaFragmentError, ModernaFragment, unknown)
        

class ModernaFragment3Tests(TestCase):
    def setUp(self):
        """
        fragment attached at 5' end of model.
        """
        self.m = RnaModel(data_type='file',data=MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
        self.s1 = ModernaStructure('file', SMALL_FRAGMENT, seq=Sequence("GCGG"))
        self.s2 = ModernaStructure('file',SMALL_FRAGMENT, seq=Sequence("GCGG"))

    def test_attributes(self):
        """Attributes of fragment are accessible."""
        three = ModernaFragment3(self.s1, anchor3=self.m['1'])
        self.assertEqual(len(three.anchor_residues), 1)
        self.assertEqual(len(three.nonanchor_residues), 3)
        self.assertTrue(str(three))

    def test_add_on_5_end(self):
        """Model sequence should change accordingly."""
        three = ModernaFragment3(self.s1, anchor3=self.m['1'])
        self.m.insert_fragment(three)
        self.assertEqual(self.m.get_sequence(), Sequence("GCGGCGGAUUUALCUCAG"))
        self.assertTrue(three.rmsd<=1.00)
        
    def test_add_on_5_end_exact(self):
        """Model sequence should change accordingly."""
        # exactly fitting example replacing original residues.
        three_overwrite = ModernaFragment3(self.s2, anchor3=self.m['6'])
        self.m.insert_fragment(three_overwrite)
        self.assertEqual(self.m.get_sequence(),Sequence("GCGUUUALCUCAG"))
        self.assertTrue(three_overwrite.rmsd<=1.00)

    def test_add_discontinuous(self):
        """Tries to add a helix to an end."""
        helix = ModernaStructure('file','test_data/rna_structures/helix.pdb')
        helix_frag = ModernaFragment3(helix, anchor3=self.m['1'], strict=False)
        self.m.insert_fragment(helix_frag)
        self.assertEqual(self.m.get_sequence(), Sequence('CCGACCUUCGGCC_GGUGGCCGAAGGGCGGAUUUALCUCAG'))


class ModernaFragment5Tests(TestCase):
    """fragment attached at 3' end of model"""
    def setUp(self):
        self.m = RnaModel(data_type='file',data=MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
        self.s1 = ModernaStructure('file', SMALL_FRAGMENT, seq=Sequence("GCGG"))
        self.s2 = ModernaStructure('file',SMALL_FRAGMENT, seq=Sequence("GCGG"))

    def test_attributes(self):
        five = ModernaFragment5(self.s1, anchor5=self.m['15'])
        self.assertEqual(len(five.anchor_residues), 1)
        self.assertEqual(len(five.nonanchor_residues), 3)
        self.assertTrue(str(five))

    def test_add_on_3_end(self):
        """Model sequence should change accordingly."""
        five = ModernaFragment5(self.s1, anchor5=self.m['15'])
        self.m.insert_fragment(five)
        self.assertEqual(self.m.get_sequence(),Sequence("GCGGAUUUALCUCAGCGG"))
        self.assertTrue(five.rmsd<=1.00)

    def test_add_on_3_end_replace(self):
        """Model sequence should change accordingly."""
        five_overwrite = ModernaFragment5(self.s2, anchor5=self.m['9'])
        self.m.insert_fragment(five_overwrite)
        self.assertEqual(self.m.get_sequence(),Sequence("GCGGAUUUACGG"))
        self.assertTrue(five_overwrite.rmsd<=1.00,0)


class ModernaFragment53Tests(TestCase):
    """
    Checks basic functionality of the ModernaFragment class.
    """
    def setUp(self):
        self.m = RnaModel(data_type='file',data=MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
        self.s1 = ModernaStructure('file', SMALL_FRAGMENT, seq=Sequence("GCGG"))
        self.s2 = ModernaStructure('file',SMALL_FRAGMENT, seq=Sequence("GCGG"))
        self.s3 = ModernaStructure('file', MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
    
    def test_attributes(self):
        """object attributes are set up correctly."""
        f1 = ModernaFragment53(self.s1, anchor5=self.m['10'], anchor3=self.m['13'])
        self.assertEqual(len(f1.struc), 4)
        self.assertEqual(f1.anchor5.fixed_resi, self.m['10'])
        self.assertEqual(f1.anchor3.fixed_resi, self.m['13'])
        self.assertEqual(f1.struc.get_sequence(),Sequence("GCGG"))
        self.assertTrue(str(f1))
        # anchor residues
        self.assertEqual(len(f1.anchor_residues), 2)
        self.assertEqual(f1.anchor_residues[0].mobile_resi.identifier, '1')
        self.assertEqual(f1.anchor_residues[1].mobile_resi.identifier, '4')
        f1.prepare_anchor_residues()
        self.assertEqual(f1.anchor_residues[0].mobile_resi.identifier, '1')
        self.assertEqual(f1.anchor_residues[1].mobile_resi.identifier, '4')
        # non-anchors
        self.assertEqual(len(f1.nonanchor_residues), 2)
        # second example
        f2 = ModernaFragment53(self.s2,anchor5=self.m['1'], anchor3=self.m['4'], new_sequence=Sequence('GL'))        
        self.assertEqual(len(f2.struc), 4)
        self.assertEqual(f2.struc.get_sequence(),Sequence("GCGG"))
        self.assertEqual(f2.new_sequence,Sequence("GL"))

    def test_get_resi_to_remove(self):
        """Should return resi identifiers between anchors and of anchors itself."""
        f1 = ModernaFragment53(self.s1, anchor5=self.m['10'], anchor3=self.m['13'])
        result = f1.get_resi_to_remove(self.m)
        self.assertEqual(result, ['10', '11', '12', '13'])

    def test_renumber(self):
        """New numbering should start at anchor5, and then go through alphabet."""
        f1 = ModernaFragment53(self.s1, anchor5=self.m['10'], anchor3=self.m['13'])
        f1.prepare_anchor_residues()
        f1.renumber()
        self.assertEqual([r.identifier for r in self.s1],['10','11','12','13'])
        # second example
        f2 = ModernaFragment53(self.s2,anchor5=self.m['1'], anchor3=self.m['4'], new_sequence=Sequence('GL'))        
        f2.prepare_anchor_residues()
        f2.renumber()
        self.assertEqual([r.identifier for r in self.s2], ['1','2','3','4'])

    def test_superimpose_fragment(self):
        """Should apply superimposition and return RMSD"""
        f1 = ModernaFragment53(self.s1, anchor5=self.m['10'],anchor3=self.m['13'])
        rmsd = f1.superimpose()
        self.assertAlmostEqual(rmsd,1.0,0) # MM: ? 1.1253267913922658
        # second fragment should fit perfectly  
        f2 = ModernaFragment53(self.s2,anchor5=self.m['1'],anchor3=self.m['4'], new_sequence=Sequence('GL'))        
        rmsd = f2.superimpose()
        self.assertAlmostEqual(rmsd,0.00)

    def test_refine(self):
        """Should change both numbers and sequence, and superimpose."""
        f2 = ModernaFragment53(self.s2,anchor5=self.m['1'],anchor3=self.m['4'], new_sequence=Sequence('GL'))        
        f2.superimpose()
        f2.prepare_anchor_residues()
        f2.renumber()
        f2.apply_seq()
        numbers = [r.identifier for r in self.s2]
        self.assertEqual(numbers,['1','1A','1B','4'])
        self.assertEqual(self.s2.get_sequence(), Sequence('GGLG'))
        
    def test_clash(self):
        """Should not clash with a given piece of structure."""
        f1 = ModernaFragment53(self.s1, anchor5=self.m['10'],anchor3=self.m['13'])
        self.assertFalse(f1.has_clashes(self.m['8':'10']))
        f2 = ModernaFragment53(self.s2,anchor5=self.m['1'],anchor3=self.m['4'], new_sequence=Sequence('GL'))        
        self.assertTrue(f2.has_clashes(self.m['1':'4']))

    def test_create_long_fragment(self):
        """Fragments with 26+ residues should be OK."""
        t = load_template(RNA_1C0A, 'B')
        resis = t['5':'45']
        long_seq = Sequence("AUAUAUAUAUGCGCGCGCGCAUAUAUAUAUGCGCGCGCG")
        struc = ModernaStructure('residues', resis)
        frag = ModernaFragment53(struc, anchor5=t['5'],anchor3=t['45'], new_sequence=long_seq)
        frag.superimpose()
        frag.prepare_anchor_residues()
        frag.renumber()
        frag.apply_seq()
        self.assertEqual(struc.get_sequence(), Sequence("CAUAUAUAUAUGCGCGCGCGCAUAUAUAUAUGCGCGCGCGG"))
    
    def test_add_numbering_letters(self):
        """Model numbers should change accordingly."""
        middle1 = ModernaFragment53(self.s1,anchor5=self.m['9'],anchor3=self.m['14'], keep=keep_nothing)
        middle1.prepare_anchor_residues()
        middle1.renumber(self.m)
        self.m.insert_fragment(middle1)
        numbers = [r.identifier for r in self.m]
        expected = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '9A', '9B','14','15']
        self.assertEqual(numbers, expected)

    def test_add_keep_numbers(self):
        """Model numbers should change accordingly."""
        middle1 = ModernaFragment53(self.s1,anchor5=self.m['9'],anchor3=self.m['14'], keep=keep_first_last)
        middle1.prepare_anchor_residues()
        middle1.renumber(self.m)
        self.m.insert_fragment(middle1)
        numbers = [r.identifier for r in self.m]
        expected = ['1', '2', '3', '4', '5', '6', '7', '8','9', '10', '13', '14','15']
        self.assertEqual(numbers, expected)

    def test_add_continuous_exact(self):
        """Model sequence should change accordingly."""
        f2 = ModernaFragment53(self.s2,anchor5=self.m['1'], anchor3=self.m['4'], new_sequence=Sequence('GL'))        
        self.m.insert_fragment(f2)
        self.assertEqual(self.m.get_sequence(), Sequence("GGLGAUUUALCUCAG"))
        self.assertAlmostEqual(f2.rmsd,0.000,2)



if __name__ == '__main__':
    main()
    
        
