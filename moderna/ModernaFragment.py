#!/usr/bin/env python
#
# ModernaFragment.py
#
# Class that supports inserting missing fragments into a model. 
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


from ModernaSuperimposer import ModernaSuperimposer
from ModernaResidue import ModernaResidue
from analyze.ClashRecognizer import ClashRecognizer
from Renumerator import LetterRenumerator, NumberRenumerator
from Errors import ModernaFragmentError
from Constants import LIR_SUPERPOSITION5, LIR_SUPERPOSITION3

#
# anchor residue classes
#
class AnchorBuildRule(object):
    """
    Instructions which atoms to copy to an anchor.
    """
    def __init__(self, base, fixed_at, mobile_at):
        """
        base - defines from which residue (fragment/model) to copy atoms.
                    It can be either 'fragment' or 'model'
        fixed_at - defines extra atoms to copy from model (atom names list)
        mobile_at - defines extra atoms to copy from fragment (atom names list)
        """
        self.base = base
        self.fixed_at = fixed_at
        self.mobile_at = mobile_at
        if self.base not in VALID_BASE:
            raise ModernaFragmentError("anchor copy option 'base'\
 must be from %s"%str(VALID_BASE))
        
VALID_BASE = ('fragment', 'model')
ALL_FROM_MODEL = AnchorBuildRule('model', [], [])
ALL_FROM_FRAGMENT = AnchorBuildRule('fragment', [], [])

        
class AnchorResidue(object):
    """
    A pair of residues 
    that can be superimposed and merged 
    for inserting fragments.
    """
    def __init__(self, fixed_resi, mobile_resi, superpos_atoms, build_rule):
        self.fixed_resi = fixed_resi
        self.mobile_resi = mobile_resi
        self.sup_atoms = superpos_atoms
        self.build_rule = build_rule
        self._validate()
        self.fixed_id = fixed_resi.identifier
        self.mobile_id = mobile_resi.identifier
        
    def __repr__(self):
        return '<Anchor residue pair: fixed %s; mobile %s>'\
            % (self.fixed_id, self.mobile_id)

    def _validate(self):
        """Checks whether data for the anchor is complete."""
        if not self.fixed_resi: 
            raise ModernaFragmentError('Expects a fixed residue as anchor.')
        if not self.mobile_resi: 
            raise ModernaFragmentError('Expects a mobile residue as anchor.')

    @property
    def fixed_superposition_atoms(self):
        """Returns atoms for superposition from the model."""
        return self.fixed_resi.get_atoms_by_names(self.sup_atoms, strict=True)

    @property
    def mobile_superposition_atoms(self):
        """Returns atoms for superposition from the fragment."""
        return self.mobile_resi.get_atoms_by_names(self.sup_atoms, strict=True)

    def _build_residue_from_two(self, resi, second_resi, atom_names):
        """Combines atoms from two residues into one."""
        prep = ModernaResidue(resi, alphabet_entry=resi.alphabet_entry)
        for atom in atom_names:
            if prep.has_id(atom):
                prep.detach_child(atom)
            prep.add(second_resi[atom])
        return prep
        
    def get_anchored_residue(self):
        """Returns a ModernaResidue that can be inserted as a anchor."""
        if self.build_rule.base == 'model': # fixed
            prep = self._build_residue_from_two(self.fixed_resi, \
                            self.mobile_resi, self.build_rule.mobile_at)
            prep.mutate(self.fixed_resi.long_abbrev) 
        elif self.build_rule.base == 'fragment': # mobile
            prep = self._build_residue_from_two(self.mobile_resi, \
                            self.fixed_resi, self.build_rule.fixed_at)
        return prep



#
# Rules for keeping residue numbers from model.
# (uses the strategy pattern)
#
def keep_nothing(frag, new_ids, struc):
    """Leave residue ids unchanged."""
    return new_ids
    
def keep_first(frag, new_ids, struc):
    """Copy first resi id after anchor5."""
    kept = struc.find_residues_in_range(frag.anchor5.fixed_id)
    if len(kept)>0:
        new_ids[1] = kept[0]
    return new_ids

def keep_last(frag, new_ids, struc):
    """Copy last resi id before anchor3."""
    kept = struc.find_residues_in_range(None, frag.anchor3.fixed_id)
    if len(kept)>0:
        new_ids[-1] = kept[-1]
    return new_ids

def keep_first_last(frag, new_ids, struc):
    """Copy two resi ids adjacent to 5'+3' anchors."""
    kept = struc.find_residues_in_range(frag.anchor5.fixed_id, \
                                        frag.anchor3.fixed_id)
    if len(kept)>1:
        new_ids[1] = kept[0]
        new_ids[-2] = kept[-1]
    return new_ids


#
# Fragment classes
#
class ModernaFragment(object):
    """
    A piece of structure that can be attached to a model.
    """
    def __init__(self, struc, new_sequence=None, \
                keep=keep_nothing, strict=True):
        self.struc = struc
        self.new_sequence = new_sequence
        self.sup = ModernaSuperimposer()
        self.rmsd = None
        self._keep_func = keep
        if strict: 
            self.check()
        self._prepared_anchors = []
        
    @property
    def anchor_residues(self):
        """Return all anchor residues."""
        return []

    @property
    def nonanchor_residues(self): 
        """Returns all non-anchor residues."""
        return list(self.struc)

    def __str__(self):
        fr_str = '%s (length %i)\n' % (self.__class__.__name__, len(self.struc))
        fr_str += '  anchor residues: %s\n' % str(self.anchor_residues)
        fr_str += '  sequence       : %s\n' % str(self.struc.get_sequence())
        fr_str += '  sec structure  : %s\n' % str(self.struc.get_secstruc())
        return fr_str
 
    def check(self):
        """Checks if the sequence given is OK."""
        if not self.struc.is_chain_continuous(): 
            raise ModernaFragmentError('Cannot create a ModernaFragment \
instance. The backbone is not continuous')
        if '.' in self.struc.get_sequence().seq_with_modifications: 
            raise ModernaFragmentError('Cannot create a ModernaFragment. \
Unknown residue is present in the chain.')
 
    def _superimpose_anchors(self, anchors, atoms):
        """Helper function that allows to select some anchors."""
        at_model = []
        at_frag = []
        for anchor in anchors:
            at_model += anchor.fixed_superposition_atoms
            at_frag += anchor.mobile_superposition_atoms
        self.sup.fixed = at_model
        self.sup.moved = at_frag
        self.sup.moved_atoms = atoms
        rmsd = self.sup.superimpose()
        return rmsd
                
    def superimpose(self): 
        """Superimpose the fragment using all anchor residues. Returns RMSD."""
        atoms = self.struc.get_all_atoms()
        self.rmsd = self._superimpose_anchors(self.anchor_residues, atoms)
        return self.rmsd
        
    def _get_resis_to_apply_seq(self):
        """Returns residues whose sequence should change."""
        return self.nonanchor_residues

    def apply_seq(self):
        """Changes seq of fr residues between anchors (NOT anchors)"""
        if self.new_sequence:
            resis = self._get_resis_to_apply_seq()
            if len(resis) != len(self.new_sequence):  
                raise ModernaFragmentError('Length of sequence (%i) does not \
match number of fragment residues to change (%i).'\
                %(len(self.new_sequence), len(resis))) 
            for resi, letter in zip(resis, self.new_sequence):
                resi.mutate(letter.long_abbrev)

    def get_resi_to_remove(self, struc):
        """
        returns residue ids from struct to be replaced 
        by the fragment plus anchors.
        """
        return []

    def fix_backbone(self): 
        """Repairs backbone breaks at all anchors."""
        pass

    def get_original_anchors(self):
        """Returns a list of unmodified anchors from fragment."""
        return [anchor.mobile_resi for anchor in self.anchor_residues]

    def prepare_anchor_residues(self):
        """
        takes care about right atoms in anchor, 
        applys model sequence for anchors, 
        does NOT change anchor numeration
        """
        self._prepared_anchors = []
        self.superimpose()
        for anchor in self.anchor_residues:
            self._prepared_anchors.append(anchor.get_anchored_residue())
      
    def _get_resis_to_renumber(self):
        """Returns a list of residues to be renumbered."""
        return self._prepared_anchors + self.nonanchor_residues

    def _get_numbers_to_renumber(self, struc):
        """Returns a list of kept ids or Nones for renumbering."""
        new_ids = [a.identifier for a in self._prepared_anchors]
        new_ids += [None] * len(self.nonanchor_residues)
        return new_ids

    def _get_prepared_gaps(self):
        """Returns a list of anchors and gap lengths required for renumbering."""
        nums = self._get_numbers_to_renumber(None)
        gaps = []
        last_n = 0
        anchor = 0
        gap_length = 0
        shift = 0
        for n in nums:
            if n is not None:
                if last_n is None:
                    gaps.append((str(int(anchor) + shift), gap_length))
                    shift += gap_length
                    gap_length = 0
                else:
                    anchor = n
            else:
                gap_length += 1
            last_n = n
        return gaps

    def _renumber_fixed_anchor_ids(self, model):
        """Changes the numbers of fixed anchors to conform the new numbering."""
        pass

    def renumber(self, struc=None, renumerator=NumberRenumerator):
        """
        changes anchor and all numeration of all residues, 
        takes care about numbers that should be kept from model.
        """
        new_resi = self._get_resis_to_renumber()
        new_ids = self._get_numbers_to_renumber(struc)
        new_ids = self._keep_func(self, new_ids, struc)
        self.struc.remove_all_residues()
        # renumber
        new_ids = renumerator(new_ids).get_identifiers()
        for new_id, resi in zip(new_ids, new_resi):
            self.struc.add_residue(resi, number=new_id)
        
    def has_clashes(self, resi_list):
        """
        Finds out whether the fragment clashes with the given list of residues.
        """
        clashrec = ClashRecognizer()
        clashes = clashrec.find_clashes_in_residues(self.nonanchor_residues+resi_list)
        return clashes
            
  
class ModernaFragment5(ModernaFragment):
    """
    Fragment connected to the model by one residue on its 5' end.
    """
    def __init__(self, struc, anchor5=None, new_sequence=None, \
                sup_r5=LIR_SUPERPOSITION5, build_rule=ALL_FROM_MODEL, \
                keep=keep_nothing, strict=True):
        ModernaFragment.__init__(self, struc, new_sequence, keep, strict)
        self.anchor5 = AnchorResidue(anchor5, list(self.struc)[0], sup_r5, build_rule)

    @property
    def anchor_residues(self):
        """Return anchors at 5' end of the fragment."""
        return [self.anchor5]
        
    @property
    def nonanchor_residues(self):
        """Returns all non-anchor residues."""
        return [r for r in self.struc][1:]
   
    def get_resi_to_remove(self, struc):
        """
        returns residue ids from struct to be replaced 
        by the fragment plus anchors.
        """
        return [self.anchor5.fixed_id] \
        + struc.find_residues_in_range(self.anchor5.fixed_id)
    
    def fix_backbone(self): 
        """Repairs backbone breaks at all anchors."""
        self.struc.fix_backbone_after_resi(self.struc[self.anchor5.fixed_id])


class ModernaFragment3(ModernaFragment):
    """
    Fragment connected to the model by one residue on its 3' end.
    """
    def __init__(self, struc, anchor3=None, new_sequence=None, \
                sup_r3=LIR_SUPERPOSITION3, build_rule=ALL_FROM_MODEL, \
                keep=keep_nothing, strict=True):        
        ModernaFragment.__init__(self, struc, new_sequence, keep, strict)
        self.anchor3 = AnchorResidue(anchor3, list(self.struc)[-1], \
                                                   sup_r3, ALL_FROM_MODEL)

    @property
    def anchor_residues(self):
        """Return anchors at 3' end of the fragment."""
        return [self.anchor3]

    @property
    def nonanchor_residues(self):
        """Returns all non-anchor residues."""
        return [r for r in self.struc][:-1]    

    def get_resi_to_remove(self, struc):
        """returns residue ids from struct to be replaced by the fragment plus anchors."""
        resis = struc.find_residues_in_range(None, self.anchor3.fixed_id)
        return resis + [self.anchor3.fixed_id]
        
    def _get_resis_to_renumber(self):
        """Returns a list of residues to be renumbered."""
        return self.nonanchor_residues + self._prepared_anchors 

    def _get_numbers_to_renumber(self, struc):
        """Returns a list of kept ids or Nones for renumbering."""
        return [None] * len(self.nonanchor_residues) + [self.anchor3.fixed_id]

    def _renumber_fixed_anchor_ids(self, model):
        gaps = self._prepared_gaps
        new_num = str(int(gaps[0][0]) + gaps[0][1] + 1)
        self.anchor3.fixed_resi = model[new_num]
        self.anchor3.fixed_id = new_num

    def fix_backbone(self): 
        """Repairs backbone breaks at all anchors."""
        self.struc.fix_backbone_before_resi(self.struc[self.anchor3.fixed_id])


class ModernaFragment53(ModernaFragment):
    """
    Fragment connected to the model by one residue on both its ends.
    """
    def __init__(self, struc, anchor5=None, anchor3=None, new_sequence=None, \
                sup_r5=LIR_SUPERPOSITION5, sup_r3=LIR_SUPERPOSITION3, \
                build_rule5=ALL_FROM_MODEL, build_rule3=ALL_FROM_MODEL, \
                keep=keep_nothing, strict=True):
        ModernaFragment.__init__(self, struc, new_sequence, keep, strict)
        self.anchor5 = AnchorResidue(anchor5, \
                                     list(self.struc)[0], sup_r5, build_rule5)
        self.anchor3 = AnchorResidue(anchor3, \
                                     list(self.struc)[-1], sup_r3, build_rule3)

    @property
    def anchor_residues(self):
        """Return anchors at 5' and 3' ends of the fragment."""
        return [self.anchor5, self.anchor3]

    @property
    def nonanchor_residues(self):
        """Returns all non-anchor residues."""
        return [r for r in self.struc][1:-1]

    def get_resi_to_remove(self, struc):
        """
        returns residue ids from struct to be replaced 
        by the fragment plus anchors.
        """
        return [self.anchor5.fixed_id] \
            + struc.find_residues_in_range(self.anchor5.fixed_id, \
                                           self.anchor3.fixed_id)\
            + [self.anchor3.fixed_id]
    
    def _get_resis_to_renumber(self):
        """Returns a list of residues to be renumbered."""
        return [self._prepared_anchors[0]] \
        + self.nonanchor_residues \
        + [self._prepared_anchors[1]]

    def _get_numbers_to_renumber(self, struc):
        """
        Returns a list of kept ids, including anchors, 
        and None for all residues to be renumbered.
        """
        return [self._prepared_anchors[0].identifier] \
            + [None] * len(self.nonanchor_residues) \
            + [self._prepared_anchors[1].identifier]

    def _renumber_fixed_anchor_ids(self, model):
        gaps = self._prepared_gaps
        new_num = str(int(gaps[0][0]) + gaps[0][1] + 1)
        self.anchor3.fixed_resi = model[new_num]
        self.anchor3.fixed_id = new_num
    
    def fix_backbone(self): 
        """Repairs backbone breaks at all anchors."""
        self.struc.fix_backbone_after_resi(self.struc[self.anchor5.fixed_id])
        self.struc.fix_backbone_before_resi(self.struc[self.anchor3.fixed_id])

