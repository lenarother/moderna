#!/usr/bin/env python
#
# test_MODERNA_all.py
#
# runs complete set of unit tests for ModeRNA.
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
import sys
from moderna.util.LogFile import log

# RNA sequences and alignments
from test_sequence import *

# data infrastructure
from test_rna_residue import RNAResidueTests
from test_rna_chain import RNAChainTests
from test_moderna_residue import ModernaResidueTests
from test_moderna_structure import ModernaStructureTests
from test_structure_library import StructureLibraryTests
from test_check_pdb import CheckPdbTests
from test_write_pdb import WritePDBTests
from test_moderna_superimposer import SuperimposerTests


# RNA modeling
from test_add_remove_modification import RemoveModificationTests, AddModificationTests, ExchangeModificationTests
from test_copy_residue import CopyResidueTests
from test_exchange_bases import ExchangeBaseTests
from test_chi_rotation import ChiRotationTests
from test_moderna_fragment import ModernaFragmentTests, \
    ModernaFragment5Tests, ModernaFragment3Tests, ModernaFragment53Tests, AnchorResidueTests
from test_secstruc_fragment import ModernaFragment2DTests, ModernaFragmentStrandTests
from test_helix import HelixTests, HelixFragmentBuilderTests
from test_search_lir import FragmentFinderTests, FragmentCandidatesTests, WriteFragmentCandidatesTests, LirHitTests
from test_lir import LirRecordTests
from test_lir_insertions import LIRInsertionTests
from test_lirdb import LirDbTests
from test_renumerator import RenumeratorTests
from test_fragment_insertion import FragmentInserterTests
from test_analyze import *
from test_builder import *
from test_isosteric import *

from test_rnamodel import BasicRnaModelTests, RetainTemplateTests, IndelQualityTests

# toplevel functions
from test_commands import CommandTests
from test_util import ValidatorTests, StrucValidatorTests
from test_commandline import CommandlineTests

if __name__ == '__main__':
    log.write_to_stderr = False
    log.raise_exceptions = True
    log.redirect_stdout()
    main()
