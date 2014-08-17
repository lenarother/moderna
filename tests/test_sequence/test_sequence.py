#!/usr/bin/env python
#
# test_builder.py
#
# runs tests of all sequence and alignment modules
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

from unittest import main
from moderna.LogFile import log

from test_alphabet import AlphabetTests, AlphabetEntryTests
from test_alignment_position import AlignmentPositionTests
from test_rnasequence import SequenceTests
from test_rnaalignment import RNAAlignmentTests, RNAAlignmentParserTests
from test_alignment import AlignmentTests
from test_alignment_matcher import AlignmentMatcherTests


if __name__ == '__main__':
    log.write_to_stderr = False
    log.raise_exceptions = True
    log.redirect_stdout()
    main()
