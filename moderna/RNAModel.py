#!/usr/bin/env python
#
# RNAModel.py
#
# Create a preliminary structure model of RNA molecule.
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

from Bio.PDB import PDBParser
from sequence.ModernaAlphabet import alphabet
from sequence.ModernaSequence import Sequence
from ModernaStructure import ModernaStructure
from Template import Template
from ModernaFragment import ModernaFragment53, ModernaFragment5, \
    ModernaFragment3, keep_first, keep_last
from FragmentInsertion import FragmentInserter
from ModernaResidue import ModernaResidue
from ModernaSuperimposer import ModernaSuperimposer
from LIR import *
from SearchLIR import FragmentFinder,  FragmentCandidates,  LirScoringOption, LirQuery,  LirHit
from Constants import PATH_TO_LIR_STRUCTURES, B_FACTOR_COPY, SINGLE_STRAND,  NUMBER_OF_FRAGMENT_CANDIDATES,  BACKBONE_ATOMS
from Errors import ModernaStructureError, ModernaSuperimposerError, \
            ModernaResidueError, ExchangeBaseError, \
            RemoveModificationError, AddModificationError,\
            RnaModelError, ModernaError
from LogFile import log


class RnaModel(ModernaStructure):
    """
    Collects preliminary structure model of RNA molecule.

    Arguments:
    * template - structure of template as a Template object (must corespond with template sequence in aligment)
    * alignment - template-target sequence aligment as a Alignment object
    """

    def __init__(self, template=None, alignment=None, model_chain_name='A', data_type=None, data=None, seq=None):
        ModernaStructure.__init__(self, data_type=data_type, data=data, chain_name=model_chain_name, seq=seq)
        #TODO: order of arguments is inconsistent. 
        self.template = template
        if template:
            self.template.set_template_numeration()
        self.alignment = alignment        
        self.s = ModernaSuperimposer()

##################################### COPYING #################################

    def copy_one_residue(self, residue, number_in_model=None, strict=True):
        """
        Copies a given residue to a model on given position.

        Arguments:
        - residue to copy (as a ModernaResidue or PDB.Residue.Residue instance)
        - position in model (by default position in previous structure)  
        """
        temp_res = ModernaResidue(residue)
        num = number_in_model or temp_res.identifier
        temp_res.set_bfactor(B_FACTOR_COPY)
        self.add_residue(temp_res, num, strict = strict)       
        log.write_message('Residue %s: residue copied from template residue %s to model.' %(num, temp_res.identifier))


    def copy_list_of_residues(self, list_of_moderna_residues, list_of_numbers_in_model=None, strict=True):
        """
        Copies list of given residues to a model on given positions (also a list)

        Arguments:
        - list of residues
        - list of positions
        """
        if not list_of_numbers_in_model:
            list_of_numbers_in_model = []
            for resi in list_of_moderna_residues:
                list_of_numbers_in_model.append(str(resi.id[1]).strip()+resi.id[2].strip())
                #TODO: use resi.identifier instead!

        if len(list_of_moderna_residues) != len(list_of_numbers_in_model): 
            raise RnaModelError('Number of given residues is different than number of given positions.')

        for resi_num in range(len(list_of_moderna_residues)):
            #TODO: try for resi,number in zip(list_of_moderna_residues,list_of_numbers_in_model):
            self.copy_one_residue(list_of_moderna_residues[resi_num], list_of_numbers_in_model[resi_num], strict=strict)


    def copy_all_residues(self, strict=True, modifications=True):
        """
        Copies all residues identical (acording alignment) in both: a template and ia  
        """
        if self.alignment and self.template:
            for ap in self.alignment.copy:
                res = self.template.template_residues[str(ap.template_position)]
                if not modifications and res.modified: pass 
                else: self.copy_one_residue(res, strict=strict)
        else: raise RnaModelError('There is no template or/and alignmnt')


    def copy_one_residue_backbone(self, residue, number_in_model=None, strict=True):
        """
        Copies a given residue to a model on given position.

        Arguments:
        - residue to copy (as a ModernaResidue or PDB.Residue.Residue instance)
        - position in model (by default position in previous structure)  
        """
        temp_res = ModernaResidue(residue)
        num = number_in_model or str(temp_res.id[1]).strip()+temp_res.id[2].strip()
        temp_res.make_backbone_only_residue()
        self.add_residue(temp_res, num, strict = strict)        
        log.write_message('Residue %s: residues backbone atoms copied from template to model.' %num)       


    def copy_all_residue_backbone(self, strict=True):
        """
        Copies all residues identical (acording alignment) in both: a template and ia  
        """
        if self.alignment and self.template:
            for ap in self.alignment.copy_backbone:
                res = self.template.template_residues[str(ap.template_position)]
                self.copy_one_residue_backbone(res, strict=strict)
        else: raise RnaModelError('There is no template or/and alignmnt')


################################### EXCHANGING ################################

    def exchange_one_base(self, residue, new_name, number_in_model=None):
        """   
        Exchanges base into given residue and copies it into a model on give position.
        
        Arguments:
        - residue
        - new residue name (A, G, C or U)
        - number in a model (by default old residue number)  
        """
        temp_resi = ModernaResidue(residue)
        old_name = temp_resi.long_abbrev
        num = number_in_model or temp_resi.identifier       
        temp_resi.exchange_base(new_name)
        self.add_residue(temp_resi, num)        
        log.write_message('Residue %s: base exchanged (%s ---> %s), residue added to model.' %(num, old_name, new_name)) 

 
    def exchange_list_of_bases(self, list_of_moderna_residues, list_of_new_names, list_of_numbers_in_model=None):
        """
        Exchanges bases in given residues list.

        Arguments:
        - list of residues
        - list with new names for residues
        - list with numbers that indicates new positions for residues in a model (by default old residues positions) 
        """
        if not list_of_numbers_in_model: 
            list_of_numbers_in_model = []
            for resi in list_of_moderna_residues:
                list_of_numbers_in_model.append(resi.identifier)

        if len(list_of_moderna_residues) != len(list_of_new_names) or len(list_of_new_names) != len(list_of_numbers_in_model):
            raise RnaModelError('Number of given residues is different than number of given positions.')
            
        for resi_num in range(len(list_of_moderna_residues)):
            self.exchange_one_base(list_of_moderna_residues[resi_num], list_of_new_names[resi_num], list_of_numbers_in_model[resi_num])


    def exchange_all_bases(self):
        """
        Exchanges all bases acording to given alignment.
        """
        if self.alignment and self.template:
            for ap in self.alignment.exchange:
                res = self.template.template_residues[str(ap.template_position)]
                name = ap.target_letter.original_base
                self.exchange_one_base(res, name)
        else: raise RnaModelError('There is no template or/and alignmnt')


############################### REMOVING MODIFICATIONS ##############################

    #TODO: add docstrings to these methods
    def remove_one_modification(self, num):
        """
        """
        old_name = self[num].long_abbrev
        self[num].remove_modification() 
        log.write_message('Residue %s: modification removed (%s ---> %s).' %(num, old_name, self[num].long_abbrev))
    

    def remove_all_modifications(self):
        """
        """
        for resi in self:
            if resi.modified: 
                self.remove_one_modification(resi.identifier)
                #TODO: check deletion of commented code
                #old_name = resi.long_abbrev
                #resi.remove_modification()
                #log.write_message('Residue %s: modification removed (%s ---> %s).' %(resi.identifier, old_name, resi.long_abbrev))


    def remove_one_modification_copy(self, residue, number_in_model):
        """
        """
        temp_resi = ModernaResidue(residue)
        num = number_in_model or temp_resi.number
        temp_resi.remove_modification()
        self.add_residue(temp_resi, str(num))
        log.write_message('Residue %s: modification removed (%s ---> %s).' %(num, residue.long_abbrev, temp_resi.long_abbrev))

    
    def remove_all_modifications_copy(self):
        """
        Removes all unnecessary modifications from model acordong given alignment.
        Copies all this residues without modification into model.
        """
        if self.alignment and self.template:
            for ap in self.alignment.remove_modifications:
                res = self.template.template_residues[str(ap.template_position)]
                temp_resi = ModernaResidue(res)
                temp_resi.remove_modification()
            
                if temp_resi != ap.target_letter:
                    temp_resi.exchange_base(ap.target_letter.original_base)
                self.add_residue(temp_resi)
                log.write_message('Residue %s: modification removed (%s ---> %s).' %(temp_resi.identifier, res.long_abbrev, temp_resi.long_abbrev))
        else: raise RnaModelError('There is no template or/and alignmnt')


############################### ADDING MODIFICATIONS ##############################

    #TODO: add docstrings
    def add_one_modification(self, num, modification_long_abbrev):
        """
        """
        old_name = self[num].long_abbrev
        self[num].add_modification(modification_long_abbrev)
        log.write_message('Residue %s: modification added (%s ---> %s).' %(num, old_name, modification_long_abbrev))
    

    def add_one_modification_copy(self, residue, modification_long_abbrev, number_in_model):
        """
        """
        temp_res = ModernaResidue(residue)
        num = number_in_model or temp_res.identifier
        temp_res.add_modification(modification_long_abbrev)
        self.add_residue(temp_res, num, False)
        log.write_message('Residue %s: modification added (%s ---> %s).' %(num, residue.long_abbrev, modification_long_abbrev))  
      
    def add_all_modifications_copy(self):
        """
        """
        if self.alignment and self.template:
            for ap in self.alignment.add_modifications:                
                temp_resi = ModernaResidue(self.template.template_residues[str(ap.template_position)]) 
                old_name = temp_resi.long_abbrev
                temp_resi.add_modification(ap.target_letter.long_abbrev)
                self.add_residue(temp_resi) 
                log.write_message('Residue %s: modification added (%s ---> %s).' %(temp_resi.identifier, old_name, temp_resi.long_abbrev))
        else: raise RnaModelError('There is no template or/and alignmnt')
      

################################### INSERTING ##################################

    def find_fragment_candidates(self, res5, res3, sequence, candidates_number=NUMBER_OF_FRAGMENT_CANDIDATES, lir_path=PATH_TO_LIR_STRUCTURES, secstruc=None):
        # en: candidate_number
        """
        Looks for lfragment candidates for missing fragment in a structure.
        Returns list of fragment candidates.

        Arguments:
        - anchor residue on 5' end as a ModernaResidue instance (residue preceding missing fragment)
        - anchor residue on 3' end as a ModernaResidue instance (residue fallowing missing fragment)
        - missing sequence
        """
        fragment_finder = FragmentFinder(res5,  res3,  sequence, self,  candidates_number, lir_path, secstruc)
        return fragment_finder.find_fragment_candidates() # FragmentCandidates instance
    

    def write_fragment_candidates(self, candidates, directory_name='fragment_candidates',  with_anchor_residues=False,  with_model=False,  log=True):
        """
        Writes all possible fragments to one pdb file.
        """
        candidates.write_fragment_candidates(directory_name,  with_anchor_residues,  with_model)

    def insert_fragment(self,  fragment):
        """Inserts fragment into model."""
        finsert = FragmentInserter()
        finsert.insert_fragment(fragment, self)

#############################################################################################

    def insert_best_fragment(self, start, stop, sequence,  candidates_number=NUMBER_OF_FRAGMENT_CANDIDATES):
        """
        Makes a LIR search to fill a gap at a given position.
        looks for possible candidates, checks the first ~20 for clashes,
        takes the best, removes anchor residue and inserts it in the model.        

        Arguments:
        - start position as a string (residue identifier)e.g. '3' - this is an identifier of the preceding anchor residue
        - stop position as above - this is an identifier of the fallowing anchor residue
        - fragment sequence as a Sequence object
        - number of fragment candidates.
        """
        log.write_message('\nSearching fragment between residue %s and residue %s.' %(str(start),str(stop)))
        fragment_finder = FragmentFinder(self[start],  self[stop],  sequence, self,  candidates_number)
        best_fragment=fragment_finder.find_fragment()
        log.write_message('\nBest fragment:\n%s'%str(best_fragment))
        self.insert_fragment(best_fragment)
        

    def insert_lir_candidate(self, candidate):
        """
        Insert the given LirHit instance into model.    
        """
        fragment = candidate.fragment_instance
        self.insert_fragment(fragment)


    def insert_all_fragments(self):
        """
        Deals with indels in the model.
        """
        if self.alignment and self.template: # KR: is template necessary?
            for fragment in self.alignment.add_fragment:
                seq=[]
                for ap in fragment:
                    if ap.target_letter: seq.append(ap.target_letter)
                start = self.template.template_residues[str(fragment[0].template_position)].identifier
                stop = self.template.template_residues[str(fragment[-1].template_position)].identifier 
                seq = seq[1:-1] 
                self.insert_best_fragment(start, stop, Sequence(seq))
        else: raise RnaModelError('There is no template or/and alignmnt')


    def _elongate_strand(self,  all_resis):
        """
        Enables to elongate single stranded helix.
        Used while applying missing 5' and 3' ends 
        when missing part is longer than provided SINGLE_STRAND fragment.
        Arguments:
        - list of curently existing residues from which the fragment will be created
        """
        front = [ModernaResidue(resi) for resi in all_resis]
        back = [ModernaResidue(resi) for resi in all_resis]
        x=1
        for resi in front:
            resi.change_number(str(x))
            x+=1
        for resi in back:
            resi.change_number(str(x))
            x+=1
            
        self.s.get_atoms([front[-2], front[-1]], BACKBONE_ATOMS,  'fixed')
        self.s.get_atoms([back[0], back[1]], BACKBONE_ATOMS,  'moved')
        self.s.moved_atoms = [at for resi in back for at in resi]
        self.s.superimpose()
        return front + back[2:]


    def add_missing_5p(self):
        """ """
        if self.alignment and len(self.alignment.add_fragment_5p)>0:
            anchor3 = [r for r in self][0]
            ap_list = self.alignment.add_fragment_5p[0]
            seq = Sequence(''.join([ap.target_letter.short_abbrev for ap in ap_list]))
            all_resis = [r for r in Template(SINGLE_STRAND,'file',  'A')]
            while len(all_resis)-1 <len(ap_list): all_resis=self._elongate_strand(all_resis)
            frag_resis = all_resis[:len(ap_list)+1]
            struc = ModernaStructure('residues',frag_resis)
            frag = ModernaFragment3(struc, anchor3=anchor3, new_sequence=seq, keep=keep_last)
            self.insert_fragment(frag)

    def add_missing_3p(self):
        """ """    
        if self.alignment and len(self.alignment.add_fragment_3p)>0:
            anchor5 = [r for r in self][-1]
            ap_list = self.alignment.add_fragment_3p[0]
            seq = Sequence(''.join([ap.target_letter.short_abbrev for ap in ap_list]))
            all_resis = [r for r in Template(SINGLE_STRAND,'file',  'A')]
            while len(all_resis)-1 <len(ap_list): all_resis=self._elongate_strand(all_resis)
            frag_resis = all_resis[:len(ap_list)+1]
            struc = ModernaStructure('residues',frag_resis)
            frag = ModernaFragment5(struc, anchor5=anchor5, new_sequence=seq, keep=keep_first)
            self.insert_fragment(frag)

########################### SECONDARY STRUCTURE MODELING ##########################

    def extend_helix(self, start5,  start3,  helix_seq):
        """ 
        e.g. helix_seq = 'AAAA_UUUU'
        """
        helix_fr = ModernaFragmentHelix(anchor5=start5, anchor3=start3, new_sequence=helix_seq, model=self)
        self.insert_fragment(helix_fr)
    
    def add_second_strand(self):
        """
        """
        pass

################################### REFINING ##################################
    
    def refine_model(self):
        """Sets B-factor and occupancy of each atom."""
        for resi in self:
            for at in resi:
                at.set_bfactor(0.0)
                at.set_occupancy(1.0)
        self.sort_residues()
       
################################### CREATE MODEL ##################################

    def apply_alignment(self):
        """
        """
        if self.alignment and self.template:
            self.copy_all_residues()
            self.copy_all_residue_backbone()
            self.exchange_all_bases()
            self.remove_all_modifications_copy()
            self.add_all_modifications_copy()
        else: raise RnaModelError('There is no template or/and alignment')
    

    def create_model(self):
        """
        """
        if self.alignment and self.template:            
            self.template.set_template_numeration()
            self.apply_alignment()
            self.insert_all_fragments()       
            self.add_missing_5p()
            self.add_missing_3p()
            self.fix_backbone()
        else: raise RnaModelError('There is no template or/and alignment')
        
