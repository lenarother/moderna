#!/usr/bin/env python
#
# test_rnamodel.py
#
# unit tests for the RNAModel class
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
from moderna.Template import Template
from moderna.RNAModel import RnaModel
from moderna.ModernaStructure import ModernaStructure
from moderna.util.Errors import ModernaStructureError, RenumeratorError
from moderna.analyze.BaseRecognizer import BaseRecognizer
from moderna.ModernaFragment import ModernaFragment53, keep_first_last
from moderna.sequence.RNAAlignment import read_alignment
from moderna.sequence.ModernaSequence import Sequence
from moderna.modifications import exchange_base
from moderna import load_model, find_fragment, copy_some_residues, \
    insert_fragment, match_template_with_alignment, clean_structure, \
    renumber_chain, match_alignment_with_model
from test_data import *

class RnaModelTests(TestCase):

    def setUp(self):
        self.a = read_alignment(MINI_ALIGNMENT)
        self.t = Template(MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
        self.m = RnaModel(self.t, self.a)
        self.seq_before = self.t.get_sequence()
        self.br = BaseRecognizer()

    def tearDown(self):
        self.a = None
        self.t = None
        self.m = None
        self.br = None
        self.seq_before = None


class BasicRnaModelTests(RnaModelTests):

    def test_create_model(self):
        """Should create the model automatically."""
        self.m.apply_alignment()
        self.m.insert_all_fragments()
        self.m.fix_backbone()
        self.assertEqual(self.m.get_sequence(),Sequence('ACUGUGAYUA[UACCU#PG'))
        
    def test_3p_extension(self):
        a = read_alignment("""> target
GCGGAUUUALCUCAGAAAAAAAAAA
> template
GCGGAUUUALCUCAG----------
        """)
        m = RnaModel(None, a, data_type='file', data=MINI_TEMPLATE)
        m.add_missing_3p()
        self.assertEqual(m.get_sequence(), a.target_seq)

    def test_long_3p_extension(self):
        a = read_alignment("""> target
GCGGAUUUALCUCAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
> template
GCGGAUUUALCUCAG----------------------------------------
        """)
        m = RnaModel(None, a, data_type='file', data=MINI_TEMPLATE)
        m.add_missing_3p()
        self.assertEqual(m.get_sequence(), a.target_seq)


    def test_5p_extension(self):
        a = read_alignment("""> target
AAAAAAAAAAGCGGAUUUALCUCAG
> template
----------GCGGAUUUALCUCAG
        """)
        m = RnaModel(None, a, data_type='file', data=MINI_TEMPLATE)
        m.add_missing_5p()
        self.assertEqual(m.get_sequence(), a.target_seq)

    def test_long_5p_extension(self):
        a = read_alignment("""> target
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCGGAUUUALCUCAG
> template
----------------------------------------GCGGAUUUALCUCAG
        """)
        # first number == 1 should not work
        m = RnaModel(None, a, data_type='file', data=MINI_TEMPLATE)
        self.assertRaises(RenumeratorError, m.add_missing_5p)
        # first number == 100 should work
        m = RnaModel(None, a, data_type='file', data=MINI_TEMPLATE)
        renumber_chain(m, 100)
        m.add_missing_5p()
        self.assertEqual(m.get_sequence(), a.target_seq)


    def test_create_model_with_gaps(self):
        """Should create the model automatically."""
        a = read_alignment(ALIGN_1B23_1QF6)
        t = Template(RNA_1B23,'file','R')
        m = RnaModel(t, a)
        m.apply_alignment()
        m.insert_all_fragments()
        self.assertEqual(m.get_sequence().seq_with_modifications.replace('_', ''), 'GCCGAUAUAGCUCAGDDGGDAGAGCAGCGCAUUCGUEAUGCGAAG7UCGUAGGTPCGACUCCUAUUAUCGGCACCA')

    def test_gaps_in_target(self):
        """Moderna should model gaps in the target, too."""
        a = read_alignment(ALIGN_TARGET_GAP)
        m = RnaModel(self.t, a)
        m.apply_alignment()
        m.insert_all_fragments()
        self.assertEqual(m.get_sequence().seq_with_modifications.replace('_', ''),'..CUGACCU#P')
        # KR: ok we need a more reasonable example.

    def test_create_model_with_unknown(self):
        """Alignment that contains unknown bases in the target."""
        a = read_alignment(MINI_ALIGNMENT_WITH_UNK)        
        m = RnaModel(self.t,a)
        m.create_model()
        self.assertEqual(m.get_sequence(),Sequence('..CUGUQQUACCU#P'))
        
    def test_doublegap_model(self):
        """Should create a model filling two gaps"""
        a = read_alignment('''> target
GGGAUAGUUCCAGABU#A
> template
GGGA-AG--CCAGABU#A
''')
        t = Template(DOUBLEGAP,'file','A')
        m = RnaModel(t, a)
        m.apply_alignment()
        m.insert_all_fragments()
        m.fix_backbone()
        self.assertEqual(m.get_sequence(),Sequence('GGGAUAGUUCCAGABU#A'))

    def test_fill_doublegap(self):
        """Should insert fragments into correct sites."""
        # prepare model
        t = Template(DOUBLEGAP,'file','A')
        m = RnaModel()
        for r in t: m.copy_residue(r)
        # insert first fragment
        struc = ModernaStructure('file',FRAGMENT1)
        f1 = ModernaFragment53(struc, anchor5=m['20'], anchor3=m['24'], \
            new_sequence=Sequence('AUA'), keep=keep_first_last)
        m.insert_fragment(f1)
        resids = [r.identifier for r in m]
        expected_resids1 = ['18', '19', '20', '21', '21A', '23', '24', '27', '28', '29', '30', '31', '32', '33', '34', '35']
        self.assertEqual(resids, expected_resids1)
        # insert second fragment
        struc = ModernaStructure('file',FRAGMENT2)
        f2 = ModernaFragment53(struc, anchor5=m['23'], anchor3=m['28'], \
            new_sequence=Sequence('AUUA'), keep=keep_first_last)
        m.insert_fragment(f2)
        resids = [r.identifier for r in m]
        expected_resids2 = ['18', '19', '20', '21', '21A', '23', '24', '24A', '24B','27', '28', '29', '30', '31', '32', '33', '34', '35']
        self.assertEqual(resids, expected_resids2)

    def test_fill_gap_numbering(self):
        """Should insert fragment with correct numbering."""
        t = Template(DOUBLEGAP,'file','A')
        m = RnaModel()
        for r in t: m.copy_residue(r)
        struc = ModernaStructure('file',FRAGMENT1)
        f1 = ModernaFragment53(struc, anchor5=m['20'], anchor3=m['24'], new_sequence=Sequence('AUA'))
        m.insert_fragment(f1)
        resids = [r.identifier for r in m]
        expected_resids1 = ['18', '19', '20', '20A', '20B', '20C', '24', '27', '28', '29', '30', '31', '32', '33', '34', '35']
        self.assertEqual(resids, expected_resids1)


    def test_oppositegap_model(self):
        """Should create a model with close gaps in the other respective sequence"""
        a = read_alignment('''> target
GGGAGAGCRUUAG-BU#A
> template
GGGAGAGCR--AGABU#A
''')
        t = Template(OPPOSITEGAP,'file','A')
        m = RnaModel(t, a)
        m.apply_alignment()
        m.insert_all_fragments()
        self.assertEqual(m.get_sequence().seq_with_modifications.replace('_', ''),'GGGAGAGCRUUAGBU#A')

    def test_model_with_hydro_template(self):
        """If the template contains hydrogens, modifications should be added."""
        t = Template(RNA_HYDRO, 'file','B')
        a = read_alignment("""> 3tra_A.pdb Z73314.1/2358-2429
UPA
> 1qru_B.pdb X55374.1/1-72
CAA
""")
        m = RnaModel(t, a)
        m.create_model()
        self.assertEqual(m.get_sequence(), Sequence('UPA'))

    def test_model_with_5p3p_ends(self):
        a = read_alignment('''> target
CAUGCGGAYYYALCUCAGGUA
> mini_template
---GCGGAUUUALCUCAG---
''')
        m = RnaModel(self.t,a)
        m.create_model()
        self.assertEqual(m.get_sequence(), Sequence('CAUGCGGAYYYALCUCAGGUA'))
        
    def test_model_with_close_gaps(self):
        t = Template('test_data/rna_structures/2du3_excerpt.ent','file','D')
        a = read_alignment('''> 1qru_B.pdb X55374.1/1-72
GCA-UUCCG
> 2du3_D.pdb
CUUUA-CCC
''')
        m = RnaModel()
        copy_some_residues([t['944'],t['946'],t['947']],m)
        lc = find_fragment(m,'946','947',Sequence('U'))
        insert_fragment(m,lc[0])
        lc = find_fragment(m,'944','946',Sequence(''))
        insert_fragment(m,lc[0])
        return #TODO: remove 
        # generate model
        m = RnaModel(t,a)
        m.create_model()
        
    def test_model_with_alignment_adjustment(self):
        """Introduces small corrections on alignment."""
        a = read_alignment("""> target
ACUGUGAYUA[UACCU#P-G
> template with small errors.
GCG7A----U.UAGCUCA_G
        """)
        t = Template(MINI_TEMPLATE,'file')
        match_template_with_alignment(t, a)
        m = RnaModel(t, a)
        m.create_model()
        self.assertEqual(m.get_sequence(), Sequence("ACUGUGAYUA[UACCU#PG"))
                                                                
    def test_number_gap(self):
        """Builds model with numbering gap in the template."""
        a = read_alignment("""> target
CCGACCUUCGGCCACCUGACAGUCCUGUGCGGGAAACCGCACAGGACUGUCAACCAGGUAAUAUAACCACCGGGAAACGGUGGUUAUAUUACCUGGUACGCCUUGACGUGGGGGAAACCCCACGUCAAGGCGUGGUGGCCGAAGGUCGG
> template
CCGACCUUCGGCCACCUGACAGUCCUGUGCGG----CCGCACAGGACUGUCAACCAGGUAAUAUAACCACCGG----CGGUGGUUAUAUUACCUGGUACGCCUUGACGUGGGG----CCCCACGUCAAGGCGUGGUGGCCGAAGGUCGG
""")
        t = Template(JMB_TEMPLATE, 'file')
        clean_structure(t)
        m = RnaModel(t, a)
        m.create_model()
        self.assertEqual(m.get_sequence().seq_without_breaks, a.target_seq)


class IndelQualityTests(RnaModelTests):
    
    def test_insert_indel_quality_1(self):
        """Insert a fragment without strand break"""
        t = Template('test_data/gaps/mini_1h4s_T_gap2.pdb', 'file','T')
        a = read_alignment('test_data/gaps/ali_gap2.fasta')
        m = RnaModel(t, a)
        self.assertEqual(m.get_sequence().seq_with_modifications.find('_'), -1)

    def test_insert_best_fragment_quality_1(self):
        """Insert best indel without a strand break"""
        m = load_model('test_data/gaps/mini_1h4s_T_gap2.pdb','T')
        cand = m.insert_best_fragment('64','67', Sequence('GAA'))
        m.fix_backbone()
        self.assertEqual(m.get_sequence().seq_with_modifications.find('_'), -1)
        
    def test_check_fragment_candidate_fitness(self):
        """Checks whether the fitness of fragment candidates increases"""
        last_fitness = None
        m = load_model('test_data/gaps/mini_1h4s_T_gap2.pdb','T')
        cand = m.find_fragment_candidates(m['64'], m['67'], Sequence('GAA'))
        for i in range(20):
            if not last_fitness: last_fitness = cand[i].score
            self.assertTrue(last_fitness <= cand[i].score)
            last_fitness = cand[i].score

    def test_gap_optimization_example_1(self):
        t = Template('test_data/gaps/mini_1h4s_T_gap2.pdb', 'file','T')
        a = read_alignment('test_data/gaps/ali_gap1.fasta')
        m = RnaModel(t, a)
        m.create_model()
            
    def test_gap_optimization_example_2(self):
        t = Template('test_data/gaps/mini_1h4s_T_gap2.pdb', 'file','T')
        a = read_alignment('test_data/gaps/ali_gap2.fasta')
        m = RnaModel(t, a)
        m.create_model()
            
    def test_gap_optimization_example_3(self):
        m = load_model('test_data/gaps/mini_1h4s_T_gap2.pdb','T')
        cand = m.find_fragment_candidates(m['63'], m['66'], Sequence('AGA'))
        hit = cand[4]
        m.insert_fragment(hit.fragment_instance)

    def test_gap_optimization_example_4(self):
        m = load_model('test_data/gaps/mini_1h4s_T_gap2.pdb','T')
        cand = m.find_fragment_candidates(m['64'], m['67'], Sequence('GAA'))
        hit = cand[6]
        m.insert_fragment(hit.fragment_instance)


class RetainTemplateTests(RnaModelTests):
    """
    Checks whether the template does not change during modeling.
    """
    def test_create_model_retain_template(self):
        """The template sequence should not change."""
        self.m.apply_alignment()
        self.m.insert_all_fragments()
        self.assertEqual(self.t.get_sequence(),self.seq_before)
    
    def test_copy_all_retain_template(self):
        """The template sequence should not change."""
        self.m.copy_all_residues()
        self.assertEqual(self.t.get_sequence(),self.seq_before)
    
    def test_exchange_single_base_retain_template(self):
        """Exchanging standard base should not change the template"""
        self.m.copy_residue(self.t['1'], '1')
        exchange_base(self.m['1'], 'U')
        self.assertEqual(self.t.get_sequence(),self.seq_before)
        self.assertEqual(self.br.identify_resi(self.t['1']), 'G')
        
    def test_exchange_all_retain_template(self):
        """Exchanging all standard bases should not change the template"""
        self.m.exchange_all_bases()
        self.assertEqual(self.t.get_sequence(),self.seq_before)
        
    def test_load_model(self):
        """A model should be loadable from a file."""
        m = RnaModel(None,None,'A','file',MINI_TEMPLATE)
        self.assertEqual(len(m), 15)

    def test_remove_modification_retain_template(self):
        """Removing modifications should not change the template"""
        self.m.remove_one_modification_copy(self.t['10'], '10')
        self.assertEqual(self.t.get_sequence(),self.seq_before)
        self.assertEqual(self.br.identify_resi(self.t['10']), 'm2G')
        
    def test_remove_all_modifications_copy_retain_template(self):
        """Removing all modifications should not change the template"""
        self.m.remove_all_modifications_copy()
        self.assertEqual(self.t.get_sequence(),self.seq_before)
        
    def test_add_modification_retain_template(self):
        """Adding modifications should not change the template"""
        self.m.add_one_modification_copy(self.t['7'], 'm66A', '7')
        self.assertEqual(self.t.get_sequence(),self.seq_before)
        self.assertEqual(self.br.identify_resi(self.t['7']), 'U')
        
    def test_add_all_modifications_copy_retain_template(self):
        """Adding all modifications should not change the template"""
        self.m.add_all_modifications_copy()
        self.assertEqual(self.t.get_sequence(),self.seq_before)
        
    def test_insert_all_fragments_retain_template(self):
        """The template sequence should not change."""
        m2=RnaModel(template=self.t, alignment=self.a, model_chain_name='A', data_type='file', data=MINI_TEMPLATE)
        m2.insert_all_fragments()
        self.assertEqual(self.t.get_sequence(),self.seq_before)

if __name__ == '__main__':
    main()
    
