#!/usr/bin/env python
"""
Unit Tests for isosteric bp structure fails
"""


from unittest import TestCase, main, skip
from moderna.Constants import DATA_PATH
from moderna.ModernaStructure import ModernaStructure
import os


class WesthoffBasepairTests(TestCase):

    @skip("untested prototype")
    def test_read_all(self):
        """All basepair fragments should be readable and contain resi 1,2"""
        path = DATA_PATH + 'Westhof_bp/'
        for fn in os.listdir(path):
            if fn.endswith('.pdb'):
                self.assertTrue(os.access(path + fn, os.F_OK))
                m = ModernaStructure('file', path + fn, 'A')
                self.assertTrue(m['1'])
                self.assertTrue(m['2'])


if __name__ == '__main__':
    main()
