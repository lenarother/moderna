#!/usr/bin/env python
#
# test_builder.py
#
# runs tests of all builders.
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
from moderna.util.LogFile import log

from test_backbone_builder import BackboneBuilderTests
from test_phosphate_builder import PhosphateBuilderTests
from test_fccd_loop_closer import FCCDLoopCloserTests
from test_coord_builder import CoordBuilderTests
from test_chi_rotation import ChiRotationTests

if __name__ == '__main__':
    log.write_to_stderr = False
    log.raise_exceptions = True
    log.redirect_stdout()
    main()
