#!/usr/bin/env python
#
# test_commands.py
#
# unit tests for moderna interface
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
from moderna import *
from moderna.LogFile import log
import os
from moderna.sequence.ModernaSequence import Sequence
from moderna.sequence.ModernaAlignment import Alignment
from moderna.ModernaFragment import ModernaFragment
from moderna.Template import Template
from moderna.RNAModel import RnaModel
from moderna.Errors import ModernaResidueError
from Bio.PDB import PDBParser
from moderna.analyze.GeometryAnalyzer import GeometryAnalyzer
from moderna.Constants import HELIX
from test_data import *

OUTPUT = 'test_data/test_model_result.ent'

#TODO: is it possible to get the results of these tests
#      as a single big log file?

class CommandTests(TestCase):
    """
    Acceptance tests for all scripting commands in Moderna.

    The order of most function parameters follows the
        T.A.M. (template, alignment, model) paradigm:
        The first parameter of a function refers to the template
        or a part of a structure which is taken,
        the second to an alignment or residue name,
        the third to the model.
    """
    def setUp(self):
        self.a = Alignment(MINI_ALIGNMENT)
        self.a2 = Alignment(MINI_ALIGNMENT_WITH_UNK)
        self.t = Template(MINI_TEMPLATE, seq=Sequence("GCGGAUUUALCUCAG"))
        self.m = RnaModel()

    def tearDown(self):
        self.a = None
        self.a2 = None
        self.t = None
        self.m = None
        
    def count_atoms(self,filename):
        # counts atoms in a PDB file
        atoms = 0
        for l in open(filename):
            if l[:4] == 'ATOM': atoms += 1
        return atoms
            
    #----------------------------------------------------------
    def test_add_modification_in_model(self):
        """Add modification to a residue already in the model."""
        copy_single_residue(self.t['5'], self.m)
        add_modification(self.m['5'], 'm1G')
        self.assertEqual(self.m.get_sequence(),Sequence('K'))

    def test_add_modification_and_copy(self):
        """Add modification and copy to model."""
        add_modification(self.t['5'],'m1G',self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('K'))
    #----------------------------------------------------------
    def test_add_all_modifications(self):
        """Add modifications from the alignment and copy."""
        add_all_modifications(self.t, self.a, self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('[_#P'))
        # KR: should the 'Y' in the loop be skipped because it is not there?
        # KR: should the residues be copied automatically?
    #----------------------------------------------------------
    def test_add_pair_to_base(self):
        """ """
        m = load_model(MINI_TEMPLATE)
        add_pair_to_base(m, '3', '20', Sequence('C'))
        self.assertEqual(m.get_sequence(), Sequence('GCGGAUUUALCUCAG_C'))
        self.assertEqual(m.get_secstruc(), '..(............)')
        m.write_pdb_file('frag.pdb')
        
    #----------------------------------------------------------
    def test_apply_alignment(self):
        """Do everything except inserting indels."""
        apply_alignment(self.t, self.a, self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('ACUGUA[UACCU#PG'))
        
    #----------------------------------------------------------
    def test_apply_all_indels(self):
        """Do inserting indels only."""
        self.m.template = self.t
        copy_some_residues(self.t['1':'15'], self.m)
        apply_all_indels(self.a, self.m)
        fix_backbone(self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('GCGGUGAYUAUUALCUCAG'))
        
    def test_analyze_geometry(self):
        result = analyze_geometry(self.t)
        self.assertTrue(isinstance(result, GeometryAnalyzer))
    #----------------------------------------------------------    
    def test_change_sequence(self):
        """Changes the entire sequence of a structure."""
        m = load_model(MINI_TEMPLATE)
        change_sequence(m, "AGCU!#7/PYLLAAA")
        self.assertEqual(m.get_sequence(), Sequence("AGCU!#7/PYLLAAA"))

    #----------------------------------------------------------    
    def test_copy_identical_residues_default(self):
        """Copy all residues that are the same,including modifications."""
        copy_identical_residues(self.t, self.a, self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('C_G_UA_CU_G'))

    def test_copy_identical_residues_no_modif(self):
        """Copy all residues that are the same, excluding modifications."""
        copy_identical_residues(self.t, self.a, self.m, modifications=False)
        self.assertEqual(self.m.get_sequence(),Sequence('C_G_UA_CU_G'))
    #----------------------------------------------------------
    def test_copy_single_residue(self):
        """Copy residue from template."""
        copy_single_residue(self.t['5'], self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('A'))
        
    def test_copy_single_residue_number(self):
        """Copy residue, allow renumbering."""
        copy_single_residue(self.t['5'], self.m, '6')
        self.assertEqual(self.m.get_sequence(),Sequence('A'))
        self.assertEqual(self.m['6'].short_abbrev,'A')
        
    def test_copy_single_residue_strict(self):
        """Copy residue, but strict option rejects ... ."""
        copy_single_residue(self.t['5'], self.m,strict=True)
        self.assertEqual(self.m.get_sequence(),Sequence('A'))
        # KR: What should the strict option do?
        
    #----------------------------------------------------------
    def test_copy_some_residues(self):
        """Copy more than one residue from template."""
        copy_some_residues([self.t['5'],self.t['7']], self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('A_U'))

    def test_copy_some_residues_range(self):
        """Copy more than one residue, allow range operator."""
        copy_some_residues(self.t['5':'8'], self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('AUUU'))
        
    def test_copy_some_residues_into_blank(self):
        """Should copy residues into an empty model without alignment"""
        t = load_template(RNA_1C0A, 'B')
        m = RnaModel()
        copy_some_residues(t['31':'35']+t['38':'42'],m)
        self.assertEqual(len(m), 10)

    def test_copy_some_residues_number(self):
        """Copy more than one residue from template."""
        copy_some_residues([self.t['5'],self.t['7']], self.m, ['4','7A'])
        self.assertEqual(self.m.get_sequence(),Sequence('A_U'))
        self.assertEqual(self.m['4'].short_abbrev,'A')
        self.assertEqual(self.m['7A'].short_abbrev,'U')

    def test_copy_some_residues_strict(self):
        """Copy more than one residue, but strict option rejects ..."""
        copy_some_residues(self.t['5':'8'], self.m, strict=True)
        self.assertEqual(self.m.get_sequence(),Sequence('AUUU'))
    #----------------------------------------------------------
    def test_create_fragment(self):
        """Loads a fragment."""
        f = create_fragment(HELIX, anchor5=self.t['3'], anchor3=self.t['9'])
        self.assertEqual(f.struc.get_sequence(), Sequence('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA_UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU'))
    #----------------------------------------------------------
    def test_create_model_empty(self):
        """No options result in an empty model."""
        m = create_model()
        self.assertEqual(m.get_sequence(),Sequence(''))
        
    def test_create_model_ali(self):
        """Apply the alignment and insert indels automatically."""
        m = create_model(self.t, self.a)
        self.assertEqual(m.get_sequence(),Sequence('ACUGUGAYUA[UACCU#PG'))

    def test_create_model_ali2(self):
        """Apply the alignment and insert indels automatically."""
        m = create_model(self.t, self.a2)
        self.assertEqual(m.get_sequence(),Sequence('..CUGUQQUACCU#P'))

    def test_create_model_ali_improve(self):
        """Correct small errors in the alignment automatically."""
        ali = Alignment(MINI_ALIGNMENT_INACCURATE)
        m = create_model(self.t, ali)
        self.assertEqual(m.get_sequence(),Sequence('ACUGUA7UACCUAPG'))

    def test_create_model_incorrect(self):
        """Giving an improper template seq should raise an exception."""
        # should also create log message.
        a = Alignment(DNA_ALIGNMENT)
        self.assertRaises(ModernaError, create_model, self.t,a)
    #----------------------------------------------------------
    def test_clean_structure(self):
        m = load_model(NASTY_PDB, 'A')
        result = clean_structure(m)
        self.assertTrue(result)
        self.assertFalse(re.search('Water residues present', str(result)))
    #----------------------------------------------------------
    def test_delete_residue(self):
        """Removes a single residue completely."""
        copy_some_residues(self.t['1':'4'],self.m)
        delete_residue('3',self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('GC_G'))
    #----------------------------------------------------------
    def test_examine_structure(self):
        m = load_model(NASTY_PDB, 'A')
        result = examine_structure(m)
        self.assertTrue(result)
        self.assertTrue(re.search('Water residues present', str(result)))
    #----------------------------------------------------------
    def test_exchange_single_base_default(self):
        """Exchange base already in the model."""
        copy_some_residues(self.t['6':'11'],self.m)
        exchange_single_base(self.m['7'],'C')
        exchange_single_base(self.m['10'],'A')
        self.assertEqual(self.m.get_sequence(),Sequence('UCUAAC'))
        
    def test_exchange_single_base_and_copy(self):
        """Exchange base and copy to the model."""
        exchange_single_base(self.t['7'],'C',self.m)
        exchange_single_base(self.t['10'],'A',self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('C_A'))
    #----------------------------------------------------------
    def test_exchange_some_bases_default(self):
        """Exchange some bases already in the model."""
        copy_some_residues(self.t['6':'11'],self.m)
        exchange_some_bases([self.m['7'],self.m['10']],'CA')
        self.assertEqual(self.m.get_sequence(),Sequence('UCUAAC'))

    def test_exchange_some_bases_range(self):
        """Exchange some bases, apply range."""
        copy_some_residues(self.t['6':'11'],self.m)
        exchange_some_bases(self.m['7':'9'],'CAG')
        self.assertEqual(self.m.get_sequence(),Sequence('UCAGLC'))
        
    def test_exchange_some_bases_and_copy(self):
        """Exchange some bases and copy from template."""
        exchange_some_bases(self.t['7':'9'],'CGA',self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('CGA'))

    def test_exchange_some_bases_and_copy_list(self):
        """Exchange some bases and copy from template."""
        exchange_some_bases([self.t['7'],self.t['9']],'CA',self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('C_A'))
    #----------------------------------------------------------
    def test_exchange_mismatches(self):
        """Exchanges all standard base mismatches and copies."""
        exchange_mismatches(self.t, self.a, self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('A_U_UA'))
    #----------------------------------------------------------
    def test_extend_helix(self):
        m = load_model(RNA_1EHZ, 'A')
        extend_helix(m, '3', '70', Sequence("AGUC_GACU"))
        self.assertEqual(m.get_sequence(), Sequence('GCGAGUCGAUUUALCUCAGDDGGGAGAGCRCCAGABU#AAYAP?UGGAG7UC?UGUGTPCG"UCCACAGAAUUGACUCGCACCA'))
        self.assertEqual(m.get_secstruc(), '(((((((((((..((((........)))).((((.........)))).....(((((.......))))))))))))))))....')
        
    #----------------------------------------------------------
    def test_find_clashes(self):
        """Finds clashes."""
        clashes = find_clashes(self.t)
        self.assertFalse(clashes)

    def test_find_clashes_empty(self):
        """Finding clashes works with empty model."""
        clashes = find_clashes(self.m)
        self.assertFalse(clashes)

    def test_find_clashes_positive(self):
        """Finding clashes returns a list of clashing residue pairs."""
        m = load_model(CLASHING_STRUCTURE)
        clashes = find_clashes(m)
        self.assertTrue(clashes)
        # check return type
        self.assertEqual(len(clashes),1)
        self.assertEqual(clashes[0][0].identifier,'3')
        self.assertEqual(clashes[0][1].identifier,'4')

    def test_find_clashes_list(self):
        """Finding clashes works with a list as input."""
        clashes = find_clashes([self.t['3'],self.t['6']])
        self.assertFalse(clashes)

    def test_find_clashes_two_structures(self):
        """Finding clashes works with residues from two structures."""
        copy_single_residue(self.t['3'], self.m)
        #clashes = find_clashes(self.t['1':'10']+ [self.m['3']])
        clashes = find_clashes([self.t['3'],self.m['3']])
        self.assertTrue(clashes)
        # check different expressions.
        clashes = find_clashes([self.t['3'],self.m['3']])
        self.assertTrue(clashes)
        
    def test_find_clashes_two_structures_negative(self):
        """Make sure that two distinct residues do not collide."""
        copy_single_residue(self.t['3'], self.m)
        clashes = find_clashes([self.t['7'],self.m['3']])
        self.assertFalse(clashes)
    #----------------------------------------------------------
    def test_find_modifications(self):
        """Finding modifications returns a dictionary."""
        mods = find_modifications(self.t)
        self.assertTrue(mods.has_key('10'))
        self.assertEqual(mods['10'].long_abbrev,'m2G')
        # KR: Check return type with MM
        
    def test_find_modifications_empty(self):
        """Finding modifications returns a dictionary."""
        mods = find_modifications(self.m)
        self.assertEqual(mods,{})
        # KR: Check return type with MM

    def test_find_modifications_advanced(self):
        """
        Tests for the modification recognizer: take 1ehz,
        run the topology matcher on it, and compare the 
        results to a manually entered list of modifications.
        """
        # check the modifications in 1ehz
        ehz_modifications = (
            ('10','L'),('16','D'),('17','D'),('26','R'),('32','B'),
            ('34','#'),('37','Y'),('39','P'),('40','?'),('46','7'),
            ('49','?'),('54','T'),('55','P'),('58','"')
            )
        t = load_template(RNA_1EHZ, 'A')
        modifications = find_modifications(t)
        self.assertEqual(len(modifications),len(ehz_modifications))
        keys = modifications.keys()
        keys.sort()
        for i in range(len(ehz_modifications)):
            self.assertEqual(keys[i],ehz_modifications[i][0])
            self.assertEqual(modifications[keys[i]].identifier, ehz_modifications[i][0])
            self.assertEqual(modifications[keys[i]].short_abbrev, ehz_modifications[i][1])
    #----------------------------------------------------------
    def test_fix_backbone(self):
        """Fixes backbone breaks."""
        m = RnaModel(data_type='file', data=FIXABLE_BACKBONE)
        fix_backbone(m)
        self.assertEqual(m.get_sequence(), Sequence('ACUGUG'))

    def test_fix_backbone_residues(self):
        """Accepts two residue numbers as parameters."""
        m = RnaModel(data_type='file', data=FIXABLE_BACKBONE)
        fix_backbone(m, '4', '5')
        self.assertEqual(m.get_sequence(), Sequence('ACUGUG'))

    #----------------------------------------------------------
    def test_find_fragment(self):
        m = load_model(RNA_HAIRPIN, 'D')
        candidates = find_fragment(m, '30', '40', Sequence("AGCUAGCU"))
        self.assertTrue(len(candidates)>0)
        self.assertTrue(isinstance(candidates[0].fragment_instance, ModernaFragment))
        
    def test_find_fragment_secstruc(self):
        m = load_model(RNA_HAIRPIN, 'D')
        candidates = find_fragment(m, '30', '40', Sequence("AGCUAGCU"), secstruc='((....))')
        self.assertTrue(len(candidates)>0)
        frag = candidates[0].fragment_instance
        finsert = FragmentInserter()
        finsert.prepare_fragment(frag, m)
        self.assertEqual(frag.struc.get_secstruc(), '(((....)))')
    
    def test_insert_fragment_cand(self):
        m = load_model(RNA_HAIRPIN, 'D')
        candidates = find_fragment(m, '30', '40', Sequence("AGCUAGCU"))
        insert_fragment(m, candidates[0])
        self.assertEqual(m.get_sequence(), Sequence("CUGAGCUAGCUC"))
        
    def test_insert_fragment_cand_secstruc(self):
        m = load_model(RNA_HAIRPIN, 'D')
        #TODO: try this example
        #candidates = find_fragment(m, '30', '40', Sequence('ACGGCCCCGU'), 20, secstruc='(((....)))')
        candidates = find_fragment(m, '30', '40', Sequence('ACCGCCCGGU'), 20, secstruc='(((....)))')
        candidates[0].fragment_instance.struc.write_pdb_file('frag2.pdb')
        insert_fragment(m, candidates[0])
        self.assertEqual(m.get_sequence(), Sequence("CUGACCGCCCGGUC"))
        self.assertEqual(m.get_secstruc(), "..((((....))))")
    
    def test_insert_two_strand_fragment(self):
        m = load_model(RNA_HAIRPIN, 'D')
        insert_two_strand_fragment(m, '30', '40', '31', '39', '195', '219', '196', '217', BULGE_MOTIF)
        fix_backbone(m)
        self.assertEqual(m.get_sequence(), Sequence("CUGCCUQUC/CGGC"))
        self.assertEqual(m.get_secstruc(), "..((.......).)")
        
    #----------------------------------------------------------
    def test_get_base_pairs(self):
        result = get_base_pairs(self.t)
        self.assertEqual(str(result['8']), '[8 WH 14]')

    def test_get_sequence(self):
        """Returning the entire sequence of a structure object."""
        seq = get_sequence(self.t)
        self.assertEqual(seq,Sequence('GCGGAUUUALCUCAG'))
        #TODO: Add test for discontinuous and strangely numbered struc.

    def test_get_sequence_chain(self):
        """Template structures should return their own sequence."""
        t = load_template(MINI_TEMPLATE, MINI_TEMPLATE_CHAIN_NAME)
        seq = get_sequence(t)
        self.assertEqual(seq,Sequence("GCGGAUUUALCUCAG"))
    
    def test_get_secstruc(self):
        """Get secondary structure in dot-bracket format."""
        h = load_template(RNA_HAIRPIN, 'D')
        ss = get_secstruc(h)
        self.assertEqual(ss, "..((.......))")

    def test_get_stacking(self):
        result = get_stacking(self.t)
        self.assertEqual(len(result), 9)
        self.assertEqual(result[0], ('1', '2', '>>'))
        
    #----------------------------------------------------------
    def test_load_alignment(self):
        """Return alignment loaded from fasta file."""
        a = load_alignment(MINI_ALIGNMENT)
        self.assertTrue(isinstance(a, Alignment))
        self.assertEqual(a.aligned_template_seq,Sequence('GCGGA----UUUALCUCAG'))
        self.assertEqual(a.aligned_target_seq,Sequence('ACUGUGAYUA[UACCU#PG'))
    #----------------------------------------------------------
    def test_load_template(self):
        """Return template loaded from PDB file."""
        t = load_template(MINI_TEMPLATE)
        self.assertTrue(isinstance(t, Template))
        self.assertEqual(t.get_sequence(),Sequence('GCGGAUUUALCUCAG'))
    #----------------------------------------------------------
    def test_load_model(self):
        """Return model loaded from PDB file."""
        m = load_model(MINI_TEMPLATE)
        self.assertTrue(isinstance(m, RnaModel))
        self.assertEqual(m.get_sequence(),Sequence('GCGGAUUUALCUCAG'))
    #----------------------------------------------------------
    def test_match_template_with_alignment(self):
        """Check if template sequence and alignment match."""
        # positive example
        self.assertTrue(match_template_with_alignment(self.t, self.a))
        # negative example
        t = load_template(CLASHING_STRUCTURE)
        self.assertFalse(match_template_with_alignment(t, self.a))
    #----------------------------------------------------------
    def test_match_model_with_alignment(self):
        """Check if model sequence and alignment match."""
        # negative example
        m = create_model()
        self.assertFalse(match_alignment_with_model(self.a, m)) 
        # positive example
        m = create_model(self.t, self.a)
        self.assertTrue(match_alignment_with_model(self.a, m)) 
    #----------------------------------------------------------
    def test_rotate_chi(self):
        """Rotates chi angle of a single base."""
        copy_some_residues(self.t['3':'5'],self.m)
        coord_before = self.m['4']['C5'].coord
        rotate_chi(self.m['4'],90)
        coord_after = self.m['4']['C5'].coord
        self.assertNotEqual(list(coord_after), list(coord_before))
                        
    #----------------------------------------------------------
    def test_remove_all_modifications(self):
        """Removes all modifications from a structure."""
        copy_some_residues(self.t['1':'15'],self.m)
        remove_all_modifications(self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('GCGGAUUUAGCUCAG'))

    def test_remove_modifications_advanced(self):
        """
        Tests for modification removal: take 1ehz, 
        remove all modifications. Run the modification recognizer. Make 
        sure that no modifications are found, but the correct number 
        of normal bases instead.
        """
        m = load_model(MINI_TEMPLATE, MINI_TEMPLATE_CHAIN_NAME)
        remove_all_modifications(m)
        modifications = find_modifications(m)
        self.assertEqual(modifications,{})
        self.assertEqual(len(m),15)
    #----------------------------------------------------------
    def test_remove_modification_default(self):
        """Removes a modification from the model."""
        copy_some_residues(self.t['1':'15'],self.m)
        remove_modification(self.m['10'])
        self.assertEqual(self.m.get_sequence(),Sequence('GCGGAUUUAGCUCAG'))
                            
    def test_remove_modification_empty(self):
        """Raises Exception if base is not modified."""
        copy_some_residues(self.t['1':'15'],self.m)
        self.assertRaises(ModernaError, remove_modification, self.m['11'])

    def test_remove_modification_and_copy(self):
        """Removes a modification, copies to model."""
        remove_modification(self.t['10'],self.m)
        self.assertEqual(self.m.get_sequence(),Sequence('G'))
    #----------------------------------------------------------
    def test_remove_mismatching_modifications(self):
        """Removes modifications that are no matches in the alignment."""
        remove_mismatching_modifications(self.t,self.a, self.m)
        # KR: !!! maybe apply_mismatching_modifications is better because of copy.
        self.assertEqual(self.m.get_sequence(),Sequence('C'))
        # KR: could add more sophisticated example.
    #----------------------------------------------------------
    def test_renumber_chain(self):
        copy_some_residues(self.t['1':'15'],self.m)
        renumber_chain(self.m,'4')
        sumnum = sum([int(r.identifier) for r in self.m])
        self.assertEqual(sumnum, 165)
    #----------------------------------------------------------
    def test_shrink_helix(self):
        m = load_model(RNA_1EHZ, 'A')
        shrink_helix(m, '3', '70', '7', '66')
        self.assertEqual(m.get_sequence(), Sequence('GCGUUALCUCAGDDGGGAGAGCRCCAGABU#AAYAP?UGGAG7UC?UGUGTPCG"UCCACAGACGCACCA'))
        self.assertEqual(m.get_secstruc(), '((((..((((........)))).((((.........)))).....(((((.......)))))))))....')
    #----------------------------------------------------------
    def test_write_model(self):
        """Creates a PDB file."""
        copy_some_residues(self.t['1':'5'],self.m)
        if os.access(OUTPUT,os.F_OK): os.remove(OUTPUT)
        write_model(self.m, OUTPUT)
        self.assertTrue(os.access(OUTPUT,os.F_OK))
        # re-read
        new = load_model(OUTPUT)
        self.assertEqual(new.get_sequence(),Sequence('GCGGA'))
        # test writing a template
        t = load_template(MINI_TEMPLATE)
        write_model(t, OUTPUT)
        new = load_model(OUTPUT)
        self.assertEqual(new.get_sequence(),Sequence('GCGGAUUUALCUCAG'))
    #----------------------------------------------------------
    def remove_frag_cand_files(self, testout):
        """Remove files with fragment candidates."""
        if os.access(testout,os.F_OK): os.remove(testout)
        if os.access(testout[:-4],os.F_OK): 
            for fn in os.listdir(testout[:-4]): os.remove(testout[:-4]+os.sep+fn)
        else:
            os.mkdir(testout[:-4])
        
    def test_write_fragment_candidates(self):
        """Fragment candidate file should be created."""
        testout = 'test_data/test_output.pdb'
        self.remove_frag_cand_files(testout)
        copy_some_residues(self.t['1':'10'],self.m)
        candidates=find_fragment(self.m,'2','8',Sequence('AGCU#'),20)
        write_fragment_candidates(candidates, testout[:-4])
        self.assertTrue(os.access(testout[:-4],os.F_OK))
    #----------------------------------------------------------
    def test_write_fragment_candidates_contents(self):
        """Fragment candidate file should contain the right number of models."""
        testout = 'test_data/test_output.pdb'
        self.remove_frag_cand_files(testout)
        copy_some_residues(self.t['1':'10'],self.m)
        candidates=find_fragment(self.m,'2','8',Sequence('AGCU#'),7)
        write_fragment_candidates(candidates, testout[:-4])
        for i, fn in enumerate(os.listdir(testout[:-4])):
            if not fn.endswith('.pdb'): continue
            p=PDBParser()
            struc = p.get_structure('cand',testout[:-4]+os.sep+fn)
            self.assertEqual(len(struc.child_list),1)
        self.assertEqual(i, 7)
    #----------------------------------------------------------
    def test_write_secstruc(self):
        """Creates a Vienna file."""
        if os.access(OUTPUT,os.F_OK): os.remove(OUTPUT)
        write_secstruc(self.t, OUTPUT)
        self.assertTrue(os.access(OUTPUT,os.F_OK))
        r = open(OUTPUT).readlines()
        self.assertTrue(re.search("\A>", r[0]))
        self.assertTrue(re.search("\AGCGGAUUUALCUCAG\n", r[1]))
        self.assertTrue(re.search("\A\.{15}\Z", r[2]))
        os.remove(OUTPUT)
    #----------------------------------------------------------

    #
    # SOME OTHER TESTS
    #    
    def test_first_level(self):
        """only checks whether a pdb file was created."""
        t = load_template(MINI_TEMPLATE, MINI_TEMPLATE_CHAIN_NAME)
        a = load_alignment(MINI_ALIGNMENT)
        m = create_model(t, a)
        if os.access(OUTPUT,os.F_OK): 
            os.remove(OUTPUT)
        write_model(m, OUTPUT)
        self.assertTrue(os.access(OUTPUT, os.F_OK))
        os.unlink(OUTPUT)

    def test_add_custom_fragment(self):
        """Fragments loaded from a file should be inserted."""
        t = load_template(RNA_1C0A, 'B')
        m = create_model()
        copy_some_residues(t['31':'35']+t['38':'42'],m)
        f = create_fragment(FRAGMENT1, anchor5=m['35'], anchor3=m['38'], chain_name='A', sequence=None)
        m.insert_fragment(f)
        self.assertTrue(m['35A'] and m['35B'])

    #----------------------------------------------------------        
    def test_read_write(self):
        """
        Tests for file processing:take 1ehz and 2-3 other PDB files. 
        Write it to a file. Read it again. Check that the structure contains 
        the right number of atoms, and that the added residues are in 
        the right position.
        """
        m = load_model(RNA_1EHZ)
        atoms_before = self.count_atoms(RNA_1EHZ)
        if os.access(OUTPUT,os.F_OK): os.remove(OUTPUT)
        write_model(m, OUTPUT)
        self.assertTrue(os.access(OUTPUT,os.F_OK))
        atoms_after = self.count_atoms(OUTPUT)
        self.assertEqual(atoms_after,atoms_before)
        os.unlink(OUTPUT)


    def test_build_model_2bte_B(self):
        """ Tests building model for 2j00_W on template 2bte_B.pdb"""
        ali_string = """> 2j00_W.pdb CP000026.1/4316393-4316321
GC-C--C--G---G--A-------U----A---G---C---U-C--AGU--------CGGU----------A-GA-G--C-------A----G---G-G----G-A--------U----------U-----G--------A-----------A----------A------------A-------------U---C--C--C-CGU--------------------------------------------------G--UC-C-U--U---G-G--U---U-C-G-----AU-----------U-----C---C--G---A---G--------U--C---C--G--G-GCA-CCA
> 2bte_B.pdb 
GC-C--G--G---G--G-------U----G---G---C---G-G--A-AU-------GGGU----------A-GACG--C-------G----C---A-U----G-A--------C_---------A-----U--------C-----------A----------U------------G-------------U---G--C--G-CAA--------------------------------------------------G-CGU-G-C--G---G-G--U---U-C-A-----AG-----------U-----C---C--C---G---C--------C--C---C--C--G-GCA-CCA"""
        ali = Alignment(ali_string)
        t = load_template(RNA_2BTE,'B')
        m = create_model(t,ali,'B')
       
    def test_build_model_1exd_B(self):
        """Tests building model for 2tra_A on template 1exd_B.pdb"""
        ali_string = """> 2tra_A
UC-C--G--U---G--A-------U----A---G---U---U-P--AAD---------GGD---------CA-GA-A--U-------G----G---G-C----G-C--------P----------U-----G--------U-----------C----------K------------C-------------G---U--G--C-CAG--------------------------------------------------A---U-?-G--G---G-G--T   ---P-C-A-----AU-----------U-----C---C--C---C---G--------U--C---G--C--G-GAG-C--
> 1ehd_B.pdb 
-G-G--G--G---U--A-------U----C---G---C---C-A--AGC---------GGU----------A-AG-G--C-------A----C---C-G----G-A--------U----------U-----C--------U-----------G----------A------------U-------------U---C--C--G-G-A--------------------------------------------------G--GU-C-G--A---G-G--U   ---U-C-G-----AA-----------U-----C---C--U---C---G--------U--A---C--C--C-CAG-CCA"""
        ali = Alignment(ali_string)
        t = load_template(RNA_1EXD,'B')
        m = create_model(t,ali,'B')

    def test_build_model_1efw_D(self):
        """ Tests building model for 1ehz_A on template 1efw_D. Should work even when 50 fragment candidates are checked"""
        ali_string = """> 1ehz_A.pdb M14856.1/1-73
GC-G--G--A---U--U-------U----A---L---C---U-C--AGD--------DGGG----------A-GA-G--C-------R----C---C-A----G-A--------B----------U-----#--------A-----------A----------Y------------A-------------P---?--U--G-GAG--------------------------------------------------7--UC-?-U--G---U-G--T   ---P-C-G-----"U-----------C-----C---A--C---A---G--------A--A---U--U--C-GCA-CCA
> 1efw_D.pdb
GG-A--G--C---G--G-------4----A---G---U---U-C--AGD--------CGGD---------DA-GA-A--U-------A----C---C-U----G-C--------C----------U-----Q--------U-----------C----------/------------C-------------G---C--A--G-GGG--------------------------------------------------7--UC-G-C--G---G-G--018U---P-C-G-----AG-----------U-----C---C--C---G---P--------C--C---G--U--U-CC-----"""
        ali = Alignment(ali_string)
        t = load_template(RNA_1EFW,'D')
        m = create_model(t,ali,'D') 
        if os.access(OUTPUT,os.F_OK): os.remove(OUTPUT)
        write_model(m, OUTPUT)
        self.assertTrue(os.access(OUTPUT,os.F_OK))


if __name__ == '__main__':
    log.write_to_stderr = False
    log.raise_exceptions = True
    log.redirect_stdout()
    main()
    
