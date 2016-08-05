#!/usr/bin/env python
"""
Unit Tests for removing modifications
"""


from unittest import main, TestCase
from moderna.ModernaStructure import ModernaStructure
from moderna.analyze.BaseRecognizer import BaseRecognizer
from moderna.modifications import remove_modification
from moderna.util.Errors import ModernaResidueError
from moderna.tests.test_data import *


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

