#!/usr/bin/env python
#
# RNAAlignment.py
#
# Parses a pairwise alignment of RNA sequences.
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

import re
from ModernaSequence import Sequence
from Constants import ANY_RESIDUE
from Errors import AlignmentError

DEFAULT_SHRINK = True

class AlignmentPosition(object):
    """
    Characterizes a single position in an alignment by:
    - target position, letter (if gap in target - None)
    - template position, letter (if gap template- None)
    - alignment position
    """

    def __init__(self, target_position=None, target_letter=None, \
                 template_position=None, template_letter=None, \
                 alignment_position=None):
        """
        Arguments:
        - target_position
        - target_letter as AlphabetEntry
        - template_position
        - template_letter as AlphabetEntry
        - alignment_position
        """
        self.target_position = target_position
        self.target_letter = target_letter # as AlphabetEntry
        self.template_position = template_position
        self.template_letter = template_letter # as AlphabetEntry
        self.alignment_position = alignment_position
        self.target_identifier = None
        self.template_identifier = None

    def __repr__(self):
        return '%s%s-%s%s' % ( \
                str(self.target_position), self.target_letter,  \
                str(self.template_position), self.template_letter, \
                )

    def get_new_notation_str(self):
        """Returns string representation using new Modomics nomenclature."""
        if self.target_letter:
            new_targ, targ_len = self.target_letter.get_new_notation_tuple()
        else:
            new_targ, targ_len = '-', 1
        if self.template_letter:
            new_temp, temp_len = self.template_letter.get_new_notation_tuple()
        else:
            new_temp, temp_len = '-', 1
        targ = ' ' * (max(targ_len, temp_len) - targ_len) + new_targ
        temp = ' ' * (max(targ_len, temp_len) - temp_len) + new_temp
        return targ, temp


    def is_identical(self):
        """
        Returns True if template and target letters are the same and are not ANY_RESIDUE.
        In other case returns False.
        """
        if not self.target_letter and not self.template_letter: 
            return True
        elif not self.target_letter or not self.template_letter: 
            return False
        elif self.target_letter.short_abbrev ==  ANY_RESIDUE or \
            self.template_letter.short_abbrev == ANY_RESIDUE: 
            return False
        elif self.target_letter == self.template_letter: 
            return True
        return False
    
    def is_different(self):
        return not self.is_identical()

    def has_template_gap(self):
        if not self.template_letter: 
            return True
        return False
        
    def has_target_gap(self):
        if not self.target_letter: 
            return True
        return False

    def has_gap(self):
        if not self.target_letter or not self.template_letter: 
            return True
        return False
            
    def is_unidentified(self):
        if self.has_gap():
           return False
        elif self.template_letter.short_abbrev == ANY_RESIDUE \
            or self.target_letter.short_abbrev == ANY_RESIDUE: 
                return True
        return False
                
    def check_underscore(self):
        #TODO: remove after kicking underscores out.
        if self.has_target_gap() or self.has_template_gap(): 
            return True
        if self.template_letter.short_abbrev == '_' \
            or self.target_letter.short_abbrev == '_':
            return False
        return True
        
    def is_mismatch(self):
        """A & N ---> this is not a mismatch, this is unidentified"""
        if self.has_gap() or self.is_unidentified():
            return False
        elif self.template_letter != self.target_letter:
            return True
        return False

    def is_transition(self):
        """purine ---> purine or pyrimidine ---> pyrimidine"""
        if self.has_gap() or self.is_unidentified():
            return False
        elif self.template_letter.original_base != self.target_letter.original_base:
            if self.target_letter.pyrimidine and self.template_letter.pyrimidine \
                or self.target_letter.purine and self.template_letter.purine:
                return True
        return False
    
    def is_transversion(self):
        """purine --->pyrmidine"""
        if self.has_gap() or self.is_unidentified():
            return False
        elif self.template_letter.original_base != self.target_letter.original_base:
            if self.target_letter.pyrimidine and self.template_letter.purine \
                or self.target_letter.purine and self.template_letter.pyrimidine: 
                return True
        return False

    def is_modified(self):
        """U ---> modified U ..."""
        if self.has_gap() or self.is_unidentified():
            return False
        elif self.is_mismatch():
            if self.target_letter.original_base == self.template_letter.original_base: 
                return True
        return False


class RNAAlignment(object):
    """
    Represents a target-template alignment 
    as a sequence of AlignmentPosition objects.

    >>> ali = RNAAlignment("> target\nAGCU-AGCU\n> template\nAGUUUUUCU")
    >>> ali[5]
    >>> ali.target_numeration[5]
    >>> ali.template_numeration[5]
    """
    #TODO: complete doctest
    def __init__(self, data, shrink=DEFAULT_SHRINK):
        """
        Creates an Alignment from a ((name1, seq1), (name2,seq2)) tuple.
        """
        self.alignment = []
        self.template_seq_name = None
        self.target_seq_name = None
        self.target_numeration = {}
        self.template_numeration = {}
        
        self._sequences = []
        self._set_alignment(data)
        
        if shrink:
            self.remove_gapped_columns()
            self.remove_excess_gaps()

    @property        
    def template_seq(self):
        """Returns template sequence as a Sequence instance."""
        return self._sequences[1].seq_without_gaps
        
    @property
    def target_seq(self):
        """Returns target sequence as a Sequence instance."""
        return self._sequences[0].seq_without_gaps

    @property
    def aligned_template_seq(self):
        return self._sequences[1]
        
    @property
    def aligned_target_seq(self):
        return self._sequences[0]

    def _get_aligned_sequences(self):
        return self._sequences
        
    def _set_aligned_sequences(self, data):
        self._set_alignment(((self.target_seq_name, data[0]), \
                                        (self.template_seq_name, data[1])))

    aligned_sequences = property(_get_aligned_sequences, \
                                 _set_aligned_sequences)

    def _set_alignment(self, data):
        """Sets name and contents of aligned sequences."""
        self.target_seq_name = data[0][0]
        self.template_seq_name = data[1][0]
        self._set_alignment_positions(data[0][1], data[1][1])
        
    def _set_alignment_positions(self, seq1, seq2):
        """Sets alignment to a list of AlignmentPosition instances."""
        self.alignment = []
        self._sequences = (seq1, seq2)
        i_alignment = 1
        i_template = 1
        i_target = 1
        ali_length = len(seq1)
        seq1 = seq1.seq_alphabet_list
        seq2 = seq2.seq_alphabet_list

        for position in range(ali_length):
            # default values for creating AlignmentPositions
            target_letter = seq1[position]
            template_letter = seq2[position]
            target_position = i_target
            template_position = i_template
            alignment_position = i_alignment
            
            if (seq1[position].short_abbrev == '-' \
                    or seq1[position].short_abbrev == '_')  \
                and (seq2[position].short_abbrev == '-' \
                     or seq2[position].short_abbrev == '_'): 
                continue

            elif seq2[position].short_abbrev == '-':
                template_position = None
                template_letter = None
                i_target += 1

            elif seq1[position].short_abbrev == '-':
                target_position = None
                target_letter = None
                i_template += 1
            else:
                i_template += 1
                i_target += 1
                
            i_alignment += 1
            apos = AlignmentPosition(
                    target_position = target_position, 
                    target_letter = target_letter,
                    template_position = template_position, 
                    template_letter = template_letter,
                    alignment_position = alignment_position
                    )
            self.alignment.append(apos)
        self._set_numeration()

    #
    # administrative methods
    #
    def __getitem__(self, args):
        """
        Allows to get AlignmentPositions using an index
        or a list of AlignmentPosition instances.
        """
        if type(args) == int:
            if args in range(1, len(self.alignment) + 1): 
                return self.alignment[args-1]
            else: 
                raise AlignmentError("There is no such position in the alignment: %i" % args)
        elif type(args) == slice:
            counter = args.start
            returned_positions = {}
            while counter <= args.stop:
                if counter <= len(self.alignment):
                    returned_positions[counter] = self.alignment[counter-1]
                    counter += 1
                else: 
                    raise AlignmentError("The alignment is shorter than %i" % counter)
            return returned_positions
        else:
            raise AlignmentError('Bad argument type.')

    def __len__(self):
        return len(self.alignment)

    def __iter__(self):
        return self.alignment.__iter__()

    def _set_numeration(self):
        """Fill in dicts target_numeration and template_numeration."""
        for apos in self.alignment:
            if apos.target_position: 
                self.target_numeration[apos.target_position] = apos
            if apos.template_position: 
                self.template_numeration[apos.template_position] = apos

    def get_identical_positions(self):
        """
        Returns positions in alignment which are identical 
        as a list of AlignmentPosition instances.
        """
        return [ap for ap in self.alignment if ap.is_identical()]

    #
    # calculation of Alignment properties
    #
    def calculate_similarity(self):
        """Returns the sequence similarity (0.0-1.0)"""
        score = 0.0
        for apos in self:
            if apos.is_unidentified(): 
                score += 0.0
            elif apos.is_identical(): 
                score += 1.0
            elif apos.is_modified(): 
                score += 0.9
            elif apos.is_transition(): 
                score += 0.5
        return score / max(len(self.target_seq), len(self.template_seq))

    def calculate_dissimilarity(self):
        """Returns the sequence dissimilarity (0.0-1.0)"""
        return 1.0 - self.calculate_similarity()

    def calculate_identity(self):
        """Returns the sequence identity (0.0-1.0)"""
        score = 0.0
        for apos in self:
            if apos.is_unidentified(): 
                score += 0.0
            elif apos.is_identical(): 
                score += 1.0
        return score / max(len(self.target_seq), len(self.template_seq))

    #
    # output and formatting
    #
    def __str__(self):
        return "\n%s\n%s\n" % (self.aligned_target_seq, \
                               self.aligned_template_seq)

    def to_fasta(self,  new_nomenclature=False, without_modifications=False):
        """Returns a string in FASTA format."""
        if new_nomenclature:
            targ_seq = ''
            temp_seq = ''
            for apos in self:
                letters = apos.get_new_notation_str()
                targ_seq += letters[0]
                temp_seq += letters[1]
            ali = '>%s\n%s\n>%s\n%s' % (self.target_seq_name, targ_seq, \
                    self.template_seq_name, temp_seq)
        elif without_modifications:
            ali = '>%s\n%s\n>%s\n%s'% (self.target_seq_name, \
                    self.aligned_target_seq.seq_without_modifications,\
                    self.template_seq_name, \
                    self.aligned_template_seq.seq_without_modifications)
        else: 
            ali = '>%s\n%s\n>%s\n%s' % (self.target_seq_name, \
                    self.aligned_target_seq.seq_with_modifications, \
                    self.template_seq_name, \
                    self.aligned_template_seq.seq_with_modifications)
        return ali

    def write_fasta_file(self, file_name, new_nomenclature=False, \
                    without_modifications=False):
        """Writes a FASTA file."""
        outfile = open(file_name, 'w')
        ali = self.to_fasta(new_nomenclature, without_modifications)
        outfile.write(ali)


    # 
    # alignment manipulations
    # 
    def remove_positions(self, seq, to_remove, keep):
        """Removes some positions from a Sequence"""
        for pos in to_remove:
            for k in keep:
                if k >= pos and len(seq) > k + 1:
                    seq = Sequence(seq[:k] + [seq[k+1], seq[k]] + seq[k+2:])
            seq = Sequence(seq[:pos] + seq[pos+1:])
        return seq

    def remove_gapped_columns(self):
        """Removes all columns with only gaps"""
        seq1, seq2 = self._sequences
        i = 0
        while i < len(seq1):
            if seq1[i].long_abbrev == '-' and seq2[i].long_abbrev == '-':
                seq1 = self.remove_positions(seq1, [i], [])
                seq2 = self.remove_positions(seq2, [i], [])
            else:
                i += 1
        self._set_alignment_positions(seq1, seq2)

    def remove_excess_gaps(self):
        """Removes excess gaps in situations like
AAAA--GGGGGG
AAAAAA--GGGG

this becomes
AAAAGGGGGG
AAAAAAGGGG
        """
        seq1, seq2 = self._sequences
        i = 1
        while i < len(seq1): # go through all positions of alignment
            assert len(seq1) == len(seq2)
            if seq1[i].long_abbrev == '-' and seq2[i-1].long_abbrev == '-':
                # found two adjacent gaps (in seq2 before seq1)
                rem1, rem2, keep = self.explore_gap_and_remove(seq1, seq2, i)
            elif seq1[i-1].long_abbrev == '-' and seq2[i].long_abbrev == '-':
                # found two adjacent gaps (in seq1 before seq2)
                rem2, rem1, keep = self.explore_gap_and_remove(seq2, seq1, i)
            else:
                rem1 = None
            # finally remove positions
            if rem1:
                seq1 = self.remove_positions(seq1, rem1, keep)
                seq2 = self.remove_positions(seq2, rem2, keep)
            i += 1
        self._set_alignment_positions(seq1, seq2)
        #TODO: should write log message.

    def explore_gap_and_remove(self, seq1, seq2, i):
        """
        Explores two adjacent gaps to the left and right.
        Finds positions to be removed, taking into account underscores.
        """
        remove_from_1 = []
        remove_from_2 = []
        keep = []

        ileft, iright = 0, 1
        while i + ileft < len(seq1) and i - iright >= 0\
        and seq1[i+ileft].long_abbrev == '-' and seq2[i-iright].long_abbrev == '-':
            if seq1[i-iright].long_abbrev == '_':
                keep.append(i - iright)
                iright += 1
            elif seq2[i + ileft].long_abbrev == '_':
                keep.append(i + ileft)
                ileft += 1
            else:
                remove_from_1.append(i)
                remove_from_2.append(i - iright)
                ileft += 1
                iright += 1
        return remove_from_1, remove_from_2, keep


class RNAAlignmentParser(object):
    """
    Parses an alignment from FASTA format.
    """
    def get_alignment(self, fasta_string, shrink=DEFAULT_SHRINK):
        """Returns a RNAAlignment object."""
        temp_aln = []

        for hit in fasta_string.split('\n>'):
            defline = hit.split('\n', 1)[0].replace('>', '')
            sequence = hit.split('\n', 1)[1:]
            sequence = re.sub("[\r\s\t]+", '', ''.join(sequence))
            temp_aln.append((defline, Sequence(sequence)))

        if not len(temp_aln[0][1]) == len(temp_aln[1][1]):
            raise AlignmentError('Template and target sequences have different size. Could not create the alignment.')
            
        alignment = RNAAlignment(temp_aln, shrink)
        #TODO: remove this check after kicking underscores out.
        for apos in alignment:
            if not apos.check_underscore():
                raise AlignmentError('Underscore symbols must be aligned to other underscores or gaps.')
        return alignment

        
    def get_alignment_from_file(self, filename, shrink=DEFAULT_SHRINK):
        """Returns a RNAAlignment object."""
        try:
            data = open(filename).read()
        except IOError:
            raise AlignmentError('File does not exist: %s '%filename)
        return self.get_alignment(data)
        
    def get_alignment_from_vienna(self, vienna_string, shrink=DEFAULT_SHRINK):
        """Place to add Secstruc stuff."""
        pass
        
    
