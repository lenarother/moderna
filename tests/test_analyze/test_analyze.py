#!/usr/bin/env python
#
# test_analyze.py
#
# runs tests of all analyzers.
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

from test_hbond_calculator import HBondCalculatorTests
from test_geometry_analyzer import GeometryAnalyzerTests
from test_geometry_statistics import AtomDefinitionTests, GeometryExpressionTests, GeometryStatisticsTests, PDBSetGeometryStatisticsTests, GeometryResultTests
from test_base_pair_calc import BasePairCalcTests
from test_secstruc_calculator import SecStrucCalculatorTests
from test_stacking_calculator import StackingCalculatorTests
from test_clash_recognizer import FindClashesTests
from test_base_recognizer import BaseRecognizerTests
from test_topology_matcher import MolParserTests, AnnotatedMoleculeTests
from test_chain_connectivity import ChainConnectivityTests

if __name__ == '__main__':
    log.write_to_stderr = False
    log.raise_exceptions = True
    log.redirect_stdout()
    main()
