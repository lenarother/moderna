#!/usr/bin/env python
#
# SecstrucFragment
#
# Enables RNA insertion of 2D fragment.
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

from ModernaStructure import ModernaStructure
from sequence.ModernaSequence import Sequence
from FragmentInsertion import FragmentInserter
from ModernaFragment import ModernaFragment, AnchorResidue, \
    ModernaFragment53, ModernaFragment5, ModernaFragment3, \
    keep_nothing, ALL_FROM_MODEL, ALL_FROM_FRAGMENT
from util.Errors import ModernaFragmentError
from Constants import HELIX_SUPERPOSITION, SINGLE_PAIR, \
                        PAIR_PURINE_SUPERPOSITION,  PAIR_PYRIMIDINE_SUPERPOSITION, \
                        BASE_PAIR_PATH, WC_BASE_PAIRS


class ModernaFragment53Strands(ModernaFragment53):
    """
    Abstract class for inserting secondary structural elements.    
    Defines two strands in the fragment: 
        strand5 (upstream) and strand3 (downstream).
    Allows for a second insertion step.
        
    Checks whether the upper part is connected (and whether on one or two sides).
    New_sequence e.g. 'AAAA_UUUU' 

                         ______
                        |         |
                        |---------|                                     upper part from model
                        |---------|
       anchor5_upper    anchor3_upper
    
frag_anchor5_upper    frag_anchor3_upper
                        |---------|                                     example bulge fragment
        strand5    |----       strand3
                        |---------|
          frag_anchor5     frag_anchor3
    
                  anchor5    anchor3
                        |---------|                                      already existing helix from model 
                        |---------|
                        
    """
    def __init__(self, struc, anchor5=None, anchor3=None, \
                anchor5_upper=None, anchor3_upper=None, \
                frag5_upper=None, frag3_upper=None, \
                new_sequence=None, \
                build_rule5=ALL_FROM_MODEL, build_rule3=ALL_FROM_MODEL, \
                superposition_atoms=HELIX_SUPERPOSITION, \
                keep=keep_nothing, model=None):
                    
        ModernaFragment53.__init__(self, struc, \
                anchor5=anchor5, anchor3=anchor3, \
                new_sequence=new_sequence, \
                sup_r5=superposition_atoms, sup_r3=superposition_atoms, \
                strict=False)
        
        self.strand5_upper_id = anchor5_upper.identifier
        self.frag5_upper_id = frag5_upper.identifier
        self.strand3_upper_id = anchor3_upper.identifier
        self.frag3_upper_id = frag3_upper.identifier
        # manage strands
        self._strand5_length = len(self._get_strand5_ids())
        self._strand3_length = len(self._get_strand3_ids())
        
        self.sup_atoms = superposition_atoms
        self.model = model
        if model == None:
            raise ModernaFragmentError('Model parameter needs to be given to fragment.')
        
    def _get_strand5_ids(self):
        """Returns initial ids for 5' strand."""
        return [r.identifier for r in self.struc[self.anchor5.mobile_id:self.frag5_upper_id]]

    def _get_strand3_ids(self):
        """Returns initial ids for 3' strand."""
        return [r.identifier for r in self.struc[self.frag3_upper_id:self.anchor3.mobile_id]]

    @property
    def strand5(self):
        """Returns residues from frag_anchor5 to frag_anchor5_upper"""
        return list(self.struc)[:self._strand5_length]

    @property
    def strand3(self):
        """Returns residues from frag_anchor5 to frag_anchor5_upper"""
        return list(self.struc)[-self._strand3_length:]
        
    @property
    def strand5_ids(self):
        """Returns ids for 5' strand."""
        return [r.identifier for r in self.strand5]

    @property
    def strand3_ids(self):
        """Returns ids for 3' strand."""
        return [r.identifier for r in self.strand3]

    def _get_new_strand5_numbers(self):
        """Returns new numbers for 5' strand."""
        return [None] * (self._strand5_length-2) + [self.strand5_upper_id] 

    def _get_new_strand3_numbers(self):
        """Returns new numbers for 3' strand."""
        return [self.strand3_upper_id] + [None] * (self._strand3_length-2) 

    def _get_numbers_to_renumber(self, struc):
        """Returns a list where residues to be renumbered are marked as None"""
        return [self.anchor5.fixed_id] \
            + self._get_new_strand5_numbers() \
            + self._get_new_strand3_numbers() \
            + [self.anchor3.fixed_id]
    
    def get_upper_fragment(self):
        """Returns fragment for the upper part."""
        return None
        
    def _get_resis_to_apply_seq(self):
        """Returns residues whose sequence should change."""
        return self.strand5[1:] \
        + self.strand3[:-1]

    def fix_backbone(self):
        """Apply the second insertion step after the first."""
        ModernaFragment53.fix_backbone(self)
        # insert upper part fragment into itself.
        finsert = FragmentInserter()
        finsert.insert_fragment(self.get_upper_fragment(), self.struc)
    

#
# Rules for keeping residue numbers from upper parts.
#
def keep_upper5(frag, new_ids, struc):
    """Copy all ids from upper part."""
    #TODO: unused so far!
    kept = struc.find_residues_in_range(frag.anchor5.fixed_id, None)
    for i, num in enumerate(kept):
        if i < len(new_ids):
            new_ids[i] = num
    return new_ids
    
def keep_upper3(frag, new_ids, struc):
    """Copy all ids from upper part."""
    #TODO: unused so far!
    kept = struc.find_residues_in_range(None, frag.anchor3.fixed_id)
    kept = kept[-len(new_ids):]
    for i, num in enumerate(kept):
        if i < len(new_ids):
            new_ids[i] = num
    return new_ids

def keep_upper53(frag, new_ids, struc):
    """Copy all ids from upper part."""
    kept = frag.struc.find_residues_in_range(frag.anchor5.fixed_id, frag.anchor3.fixed_id)
    for i, num in enumerate(kept):
        if i < len(new_ids):
            new_ids[i+1] = num
    return new_ids


class ModernaFragment553(ModernaFragment53Strands):
    """
    Fragment with a third anchor in upper5 position.
    """
    def _get_new_strand3_numbers(self):
        """Returns new numbers for 3' strand."""
        return [None] * (len(self.strand3)-2) 

    def get_upper_fragment(self):
        """Returns fragment for the upper part."""
        struc = ModernaStructure('residues', self.model[self.strand5_upper_id:self.anchor3.fixed_id])
        frag = ModernaFragment5(struc, anchor5=self.struc[self.strand5_upper_id], \
                new_sequence=None, \
                sup_r5=self.sup_atoms, build_rule=ALL_FROM_MODEL, \
                keep=keep_upper5, strict=True)
        return frag


class ModernaFragment533(ModernaFragment53Strands):
    """
    Fragment with a third anchor in upper3 position.
    """
    def _get_new_strand5_numbers(self):
        """Returns new numbers for 5' strand."""
        return [None] * (len(self.strand5)-1)

    def get_upper_fragment(self):
        """Returns fragment for the upper part."""
        struc = ModernaStructure('residues', self.model[self.anchor5.fixed_id:self.strand3_upper_id])
        frag = ModernaFragment3(struc, anchor3=self.struc[self.strand3_upper_id], \
                new_sequence=None, \
                sup_r3=self.sup_atoms, build_rule=ALL_FROM_MODEL, \
                keep=keep_upper3, strict=True)
        return frag


class ModernaFragment2D(ModernaFragment53Strands):
    """
    Fragment with two extra anchors in upper positions.
    """
    def _get_resis_to_apply_seq(self):
        """Returns residues whose sequence should change."""
        return self.strand5[1:-1] + self.strand3[1:-1]
    
    def get_upper_fragment(self):
        """Returns fragment for the upper part."""
        struc = ModernaStructure('residues', self.model[self.strand5_upper_id:self.strand3_upper_id])
        frag = ModernaFragment53(struc, anchor5=self.struc[self.strand5_upper_id], \
                anchor3=self.struc[self.strand3_upper_id], \
                new_sequence=None, sup_r5=self.sup_atoms, \
                build_rule5=ALL_FROM_FRAGMENT, \
                build_rule3=ALL_FROM_FRAGMENT, \
                keep=keep_upper53, strict=False)
        return frag


class ModernaFragmentStrand(ModernaFragment):
    """
    Single strand one residue fragment
    """
    def __init__(self, data_file=SINGLE_PAIR, anchor=None, \
                  new_sequence=None, identifier=None, strict=False):
        struc = ModernaStructure('file', data_file)
        if new_sequence:
            base = WC_BASE_PAIRS[new_sequence[0].short_abbrev]
            list(struc)[0].mutate(base)
        ModernaFragment.__init__(self, struc, new_sequence=new_sequence, keep=keep_nothing, strict=strict)
        
        sup_atoms = self.__choose_superposition_atoms(anchor)
        self.anchor = AnchorResidue(anchor, list(self.struc)[0], \
                            superpos_atoms = sup_atoms, build_rule=ALL_FROM_MODEL)
        
        self.new_sequence = new_sequence # one letter only
        self.identifier = identifier # new identifier for missing strand

    @property
    def anchor_residues(self):
        """Return all anchor residues."""
        return [self.anchor]

    @property
    def nonanchor_residues(self): 
        """Returns all non-anchor residues."""
        return [list(self.struc)[1]]
        
    def __choose_superposition_atoms(self, resi):
        """Returns atoms for superimposing the base pair."""
        if resi.purine: 
            return PAIR_PURINE_SUPERPOSITION
        elif resi.pyrimidine: 
            return PAIR_PYRIMIDINE_SUPERPOSITION
            
    def _get_numbers_to_renumber(self, struc):
        """Returns numbers for the inserted fragment."""
        new_ids = [self.anchor.fixed_resi.identifier, self.identifier]
        return new_ids
        
    def apply_sequence(self):
        """Changes sequence of the added nucleotide."""
        if self.new_sequence:
            list(self.struc)[-1].mutate(self.new_sequence.long_abbrev)
    
    def get_resi_to_remove(self, struc):
        """Returns residue ids to be removed from the model."""
        return [self.anchor.fixed_id] \
        + struc.find_residues_in_range(self.identifier, self.identifier)


# experimental function added by KR
def construct_noncanonical_pair(model, resid1, resid2, bptype):
    """
    Builds a noncanonical base pair on the first residue.
    bptype should be e.g. 'cWH_GA'
    """
    #TODO: add tests
    if model.moderna_residues.has_key(resid2):
        model.remove_residue(resid2)
    bptype += '.pdb'
    frag = ModernaFragmentStrand(data_type='file', data=BASE_PAIR_PATH+bptype, anchor=model[resid1],  identifier=resid2)
    print [r for r in frag]
    model.insert_fragment(frag)
