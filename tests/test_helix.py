#!/usr/bin/env python
"""
Unit Tests for ModernaStructure.ModernaFragment
"""

from unittest import main, TestCase
from moderna.Helix import HelixBuilder, HelixFragmentBuilder
from moderna.sequence.ModernaSequence import Sequence
from moderna.FragmentInsertion import FragmentInserter
from moderna.util.Errors import ModernaFragmentError
from moderna.Constants import HELIX
from moderna import load_template, create_model, copy_some_residues, create_fragment, load_model
from moderna.tests.test_data import *
import re


class HelixTests(TestCase):
    """
    Tests for helix building.
    """
    def test_build_helix(self):
        """A simple helix is built."""
        hb = HelixBuilder()
        helix = hb.build_helix(Sequence('GGG_CCC'))
        self.assertEqual(helix.get_sequence(), Sequence('GGG_CCC'))
        self.assertEqual(helix.get_secstruc(), "((()))")
        self.assertEqual(len(helix.strand5), 3)
        self.assertEqual(len(helix.strand3), 3)
        self.assertEqual([r.identifier for r in helix.strand5], ['1', '2', '3'])
        self.assertEqual([r.identifier for r in helix.strand3], ['401', '402', '403'])

    def test_lengths(self):
        """Helices from length 1 to 9 bp can be created."""
        UP  = "AUCGAAUUCCGG"
        DOWN = "CCGGAAUUCGAU"
        hb = HelixBuilder()
        for bp in range(1, 9):
            seq = Sequence(UP[:bp]+'_'+DOWN[-bp:])
            ss = '('*(bp) + ')'*(bp)
            helix = hb.build_helix(seq)
            self.assertEqual(len(helix), bp*2)
            self.assertEqual(helix.get_sequence(), seq)
            self.assertEqual(helix.get_secstruc(), ss)

    def test_long_helix(self):
        """Helices with length 88 can be created."""
        seq = Sequence("AUCGAAUUCCGGAAAAAGGGGGAUCGAAUUCCGGAAAAAGGGGG_CCCCCUUUUUCCGGAAUUCGAUCCCCCUUUUUCCGGAAUUCGAU")
        ss = '('*44 + ')'*44
        hb = HelixBuilder()
        helix = hb.build_helix(seq)
        self.assertEqual(helix.get_sequence(), seq)
        self.assertEqual(helix.get_secstruc(), ss)
        

class HelixFragmentBuilderTests(TestCase):
    
    def setUp(self):
        self.rna = load_model(RNA_HAIRPIN, 'D')
        self.hfb = HelixFragmentBuilder()
    
    def test_get_fragment(self):
        """Create a helical fragment."""
        frag = self.hfb.build_fragment(anchor5=self.rna['30'], anchor3=self.rna['40'], \
                sequence=Sequence('GGG_CCC'), model=self.rna)
        self.assertEqual(len(frag.struc), 10)
        self.assertEqual(frag.struc.get_sequence(), Sequence('GGGGC_GCCCC'))
        # test anchor residues
        self.assertEqual(frag.anchor5.fixed_resi, self.rna['30'])
        self.assertEqual(frag.anchor3.fixed_resi, self.rna['40'])
        self.assertEqual(frag.strand5_upper_id, '31')
        self.assertEqual(frag.strand3_upper_id, '39')
        self.assertEqual(frag.frag5_upper_id, '5')
        self.assertEqual(frag.frag3_upper_id, '401')
        # test properties
        self.assertEqual(frag.strand5_ids, ['1', '2', '3','4', '5'])
        self.assertEqual(frag.strand3_ids, ['401', '402', '403', '404', '405'])
        
    def test_insert_blunt_fragment(self):
        """Helices can be built at blunt ends."""
        model = load_model(HELIX)
        frag = self.hfb.build_fragment(anchor5=model['40'], anchor3=model['41'], \
                sequence=Sequence('AUCGAAUU_AAUUCGAU'))
        finsert = FragmentInserter()
        finsert.insert_fragment(frag, model)
        self.assertEqual(model.get_sequence(), Sequence('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUCGAAUU_AAUUCGAUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU'))
        self.assertEqual(model.get_secstruc(), '(((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))')

    def test_preparations(self):
        """After preparing the fragment has the right sequence and secstruc."""
        frag = self.hfb.build_fragment(anchor5=self.rna['30'], anchor3=self.rna['40'], \
                                       sequence=Sequence('GGG_CCC'), model=self.rna)
        finsert = FragmentInserter()
        finsert.prepare_fragment(frag, self.rna)
        frag.fix_backbone()
        self.assertEqual(frag.struc.get_sequence(), Sequence('GGGGCCUQUC/CGCCCC'))
        self.assertEqual(frag.struc.get_secstruc(), '(((((.......)))))')
        numbers = [r.identifier for r in frag.struc]
        self.assertEqual(numbers, ['30','30A','30B','30C','31','32', '33', \
            '34','35','36','37','38','39','39A','39B','39C','40'])

    def test_insertion(self):
        """Insert helical fragment into a model."""
        frag = self.hfb.build_fragment(anchor5=self.rna['30'], anchor3=self.rna['40'], \
                                       sequence=Sequence('GGG_CCC'), model=self.rna)
        self.rna.insert_fragment(frag)
        self.assertEqual(self.rna.get_sequence(), Sequence('CUGGGGCCUQUC/CGCCCC'))
        self.assertEqual(self.rna.get_secstruc(), '..(((((.......)))))')

    def test_extend_blunt_pair(self):
        """Extends single BP to helix even if numbers overlap."""
        pair = create_model()
        copy_some_residues([self.rna['30'],self.rna['40']], pair)
        pair.renumber_residue('40', '332')
        pair.renumber_residue('30', '21')
        seq = Sequence('GGGAUUUCGAAACCCAAGGUG_UACCGAGGAUGUAGGAAUUUC')
        frag = self.hfb.build_fragment(anchor5=pair['21'], anchor3=pair['332'], \
                                       sequence=seq, model=pair)
        pair.insert_fragment(frag)
        self.assertEqual(pair.get_sequence(), Sequence("GGGGAUUUCGAAACCCAAGGUG_UACCGAGGAUGUAGGAAUUUCC"))
    
    def test_insert_long_helix(self):
        """Helices with length 40 can be inserted."""
        seq = Sequence("AUCGAAUUCCGGAAAAAGGGGG_CCCCCUUUUUCCGGAAUUCGAU")
        ss = '('*18 + ')'*18
        self.rna.renumber_residue('40', '140')
        frag = self.hfb.build_fragment(anchor5=self.rna['30'], anchor3=self.rna['140'], \
                                       sequence=seq, model=self.rna)
        self.rna.insert_fragment(frag)
        self.assertEqual(self.rna.get_sequence(), Sequence('CUGAUCGAAUUCCGGAAAAAGGGGGCCUQUC/CGCCCCCUUUUUCCGGAAUUCGAUC'))


if __name__ == '__main__':
    main()
    
