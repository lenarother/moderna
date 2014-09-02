#!/usr/bin/env python
#
# test_alphabet.py
#
# unit tests for removing modifications
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
from moderna.ModernaStructure import ModernaStructure
from moderna.analyze.BaseRecognizer import BaseRecognizer
from moderna.modifications import remove_modification
from moderna.util.Errors import ModernaResidueError
from test_data import *


class RemoveModificationTests(TestCase):
    """Makes sure modifications can be removed."""

    def test_remove(self):
        struc = ModernaStructure('file', MINI_TEMPLATE)
        r10 = struc['10']
        self.assertEqual(BaseRecognizer().identify_resi(r10), 'm2G')
        remove_modification(r10)
        self.assertEqual(BaseRecognizer().identify_resi(r10), 'G')

    def test_remove_deoxy(self):
        struc = ModernaStructure('file', DNA_WITH_MISMATCH, 'E')
        r10 = struc['10']
        self.assertEqual(BaseRecognizer().identify_resi(r10), 'dG')
        remove_modification(r10)
        self.assertEqual(BaseRecognizer().identify_resi(r10), 'G')

    def test_remove_modification_empty(self):
        """Raises an Exception if base is not modified."""
        struc = ModernaStructure('file', MINI_TEMPLATE)
        self.assertRaises(ModernaResidueError, remove_modification, struc['11'])



if __name__ == '__main__':
    main()

