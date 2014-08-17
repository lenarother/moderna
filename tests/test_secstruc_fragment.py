#!/usr/bin/env python
#
# test_moderna_fragment.py
#
# unit tests for ModernaStructure.ModernaFragment
#
__author__ = "Magdalena Musielak, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Musielak"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Prototype"

from unittest import main, TestCase
from moderna.SecstrucFragment import ModernaFragment2D, ModernaFragmentStrand
from moderna.RNAModel import RnaModel
from moderna.ModernaStructure import ModernaStructure
from moderna.sequence.ModernaSequence import Sequence
from moderna.FragmentInsertion import FragmentInserter
from moderna.Errors import ModernaFragmentError
from moderna.Constants import HELIX
from moderna import fix_backbone
from moderna import load_template, create_model, copy_some_residues, create_fragment, load_model
from test_data import *
import re


class ModernaFragment2DTests(TestCase):
    
    def setUp(self):
        self.rna = load_model(RNA_HAIRPIN, 'D')
        self.motif = ModernaStructure('file', BULGE_MOTIF)
        
    def test_init(self):
        """Initialize 2D fragment with a structure and four anchor residues."""
        mf = ModernaFragment2D(self.motif, \
                        anchor5=self.rna['30'], anchor3=self.rna['40'], \
                        anchor5_upper=self.rna['31'], anchor3_upper=self.rna['39'], \
                        frag5_upper=self.motif['196'], frag3_upper=self.motif['217'], \
                        model=self.rna)
                        
        self.assertEqual(mf.anchor5.fixed_resi, self.rna['30'])
        self.assertEqual(mf.anchor3.fixed_resi, self.rna['40'])
        self.assertEqual(mf.anchor5.mobile_resi, self.motif['195'])
        self.assertEqual(mf.anchor3.mobile_resi, self.motif['219'])
        self.assertEqual(mf.strand5_upper_id, '31')
        self.assertEqual(mf.strand3_upper_id, '39')
        self.assertEqual(mf.frag5_upper_id, '196')
        self.assertEqual(mf.frag3_upper_id, '217')
        # test properties
        self.assertEqual(mf.strand5_ids, ['195', '196'])
        self.assertEqual(mf.strand3_ids, ['217', '218', '219'])
        self.assertEqual([r.identifier for r in mf.strand5], mf.strand5_ids)
        self.assertEqual([r.identifier for r in mf.strand3], mf.strand3_ids)

    def test_secstruc(self):
        """Fragment has the right secondary structure."""
        mf = ModernaFragment2D(self.motif, \
                        anchor5=self.rna['30'], anchor3=self.rna['40'], \
                        anchor5_upper=self.rna['31'], anchor3_upper=self.rna['39'], \
                        frag5_upper=self.motif['196'], frag3_upper=self.motif['217'], \
                        model=self.rna)
        self.assertEqual(mf.struc.get_secstruc(), "(().)")

    def test_renumber(self):
        """Both strands are numbered."""
        mf = ModernaFragment2D(self.motif, \
                        anchor5=self.rna['30'], anchor3=self.rna['40'], \
                        anchor5_upper=self.rna['31'], anchor3_upper=self.rna['39'], \
                        frag5_upper=self.motif['196'], frag3_upper=self.motif['217'], \
                        model=self.rna)
        mf.prepare_anchor_residues()
        self.assertEqual([r.identifier for r in mf._get_resis_to_renumber()], ['30', '196', '217', '218', '40'])
        self.assertEqual(mf._get_numbers_to_renumber(self.rna), ['30', '31','39', None, '40'])
        mf.renumber(self.rna)
        numbers = [r.identifier for r in mf.struc]
        self.assertEqual(numbers, ['30','31','39','39A','40'])

    def test_numeration_with_model(self):
        """Fragments are renumbered on both strands with everything in between"""
        mf = ModernaFragment2D(self.motif, \
                        anchor5=self.rna['30'], anchor3=self.rna['40'], \
                        anchor5_upper=self.rna['31'], anchor3_upper=self.rna['39'], \
                        frag5_upper=self.motif['196'], frag3_upper=self.motif['217'], \
                        model=self.rna)
        mf.prepare_anchor_residues()
        mf.renumber(self.rna)
        numbers = [r.identifier for r in mf.struc]
        self.assertEqual(numbers, ['30', '31', '39', '39A', '40'])

    def test_seq_and_secstruc(self):
        """After preparing everything, the fragment has the right sequence and secstruc"""
        frag = ModernaFragment2D(self.motif, \
                        anchor5=self.rna['30'], anchor3=self.rna['40'], \
                        anchor5_upper=self.rna['31'], anchor3_upper=self.rna['39'], \
                        frag5_upper=self.motif['196'], frag3_upper=self.motif['217'], \
                        model=self.rna)
        finsert = FragmentInserter()
        finsert.prepare_fragment(frag, self.rna)
        frag.fix_backbone()
        self.assertEqual(frag.struc.get_sequence(), Sequence('GCCUQUC/CGGC'))
        self.assertEqual(frag.struc.get_secstruc(), '((.......).)')
        numbers = [r.identifier for r in frag.struc]
        self.assertEqual(numbers, ['30','31','32', '33', \
            '34','35','36','37','38','39','39A','40'])
        
    def test_eliminate_middle(self):
        """Builds a fragment to shorten the connection between two anchor pairs"""
        helix = ModernaStructure('file', HELIX, 'A')
        helix = RnaModel(data_type='residues', data=helix['1':'10']+helix['71':'80'])
        frag = ModernaFragment2D(self.motif, \
                        anchor5=helix['2'], anchor3=helix['79'], \
                        anchor5_upper=helix['8'], anchor3_upper=helix['73'], \
                        frag5_upper=self.motif['196'], frag3_upper=self.motif['217'], \
                        new_sequence=Sequence('C'), \
                        model=helix)
        finsert = FragmentInserter()
        finsert.prepare_fragment(frag, helix)
        frag.fix_backbone()
        self.assertEqual(frag.struc.get_sequence(), Sequence('AAAA_UUUCU'))

    def test_insert(self):
        """Inserts 2D fragment into a model."""
        mf = ModernaFragment2D(self.motif, \
                        anchor5=self.rna['30'], anchor3=self.rna['40'], \
                        anchor5_upper=self.rna['31'], anchor3_upper=self.rna['39'], \
                        frag5_upper=self.motif['196'], frag3_upper=self.motif['217'], \
                        new_sequence=Sequence('C'), \
                        model=self.rna)
        self.rna.insert_fragment(mf)
        self.assertEqual(self.rna.get_sequence(), Sequence('CUGCCUQUC/CGCC'))
        self.assertEqual(self.rna.get_secstruc(), '..((.......).)')
        
    def test_insert_eliminate(self):
        """Inserts a 2D fragment into a model."""
        helix = load_model(HELIX, 'A')
        helix = RnaModel(None, None, 'A',  'residues', helix['1':'8']+helix['73':'80'])
        mf = ModernaFragment2D(self.motif, \
                        anchor5=helix['7'], anchor3=helix['74'], \
                        anchor5_upper=helix['8'], anchor3_upper=helix['73'], \
                        frag5_upper=self.motif['196'], frag3_upper=self.motif['217'], \
                        new_sequence=Sequence('C'), \
                        model=helix)
        helix.insert_fragment(mf)
        self.assertEqual(helix.get_sequence(), Sequence('AAAAAAAA_UCUUUUUUU'))
        self.assertEqual(helix.get_secstruc(), '(((((((().)))))))')


class ModernaFragmentStrandTests(TestCase):

    def setUp(self):
        self.rna = load_model(MINI_TEMPLATE, 'A')
    
    def test_init(self):
        mf  = ModernaFragmentStrand(anchor=self.rna['4'], identifier='20')
        #data_type='file', data=SINGLE_PAIR,  chain_name='A', anchor=None, new_sequence=None, identifier=None,  seq=None,  strict=False,  superposition_atoms=PAIR_SUPERPOSITION):
        self.assertEqual(mf.anchor.fixed_resi.identifier, '4')
        
    def test_prepare(self):
        """Prepares a single base"""
        mf  = ModernaFragmentStrand(anchor=self.rna['4'], identifier='20')
        mf.superimpose()
        mf.prepare_anchor_residues()
        mf.renumber(self.rna)
        resis = list(mf.struc)
        self.assertEqual(resis[1].identifier, '20')
        self.assertEqual(str(resis[1].get_bp(self.rna['4'])), '20 +/+ 4')
        
    def test_insert(self):
        mf  = ModernaFragmentStrand(anchor=self.rna['4'], identifier='20', new_sequence=Sequence('C'))
        self.rna.insert_fragment(mf)
        self.assertEqual(self.rna.get_sequence(), Sequence('GCGGAUUUALCUCAG_C'))
        self.assertEqual(self.rna.get_secstruc(), '...(...........)')

    def test_insert_many(self):
        SERIES = [
                  ('2', '24', 'G'), 
                  ('3', '23', 'C'), 
                  ('4', '22', 'C'), 
                  ('5', '21', 'U'), 
                  ('6', '20', 'A'), 
                  ]
        for num, ident, seq in SERIES:
            mf  = ModernaFragmentStrand(anchor=self.rna[num], identifier=ident, new_sequence=Sequence(seq))
            self.rna.insert_fragment(mf)
        self.rna.fix_backbone()
        self.assertEqual(self.rna.get_secstruc(), '.(((((.........)))))')
        self.assertEqual(self.rna.get_sequence(), Sequence('GCGGAUUUALCUCAG_AU_CCG'))


if __name__ == '__main__':
    main()
    
        
