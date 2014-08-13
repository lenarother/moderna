#!/usr/bin/env python
#
# test_westhof_bp.py
#
# unit tests for isosteric bp structure fails
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Pawel Skiba, Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Pawel Skiba"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Prototype"


from unittest import TestCase, main
from moderna.Constants import DATA_PATH
from moderna.ModernaStructure import ModernaStructure
import os

class WesthoffBasepairTests(TestCase):
    
    def test_read_all(self):
        """All basepair fragments should be readable and contain resi 1,2"""
        path = DATA_PATH+'Westhof_bp/'
        for fn in os.listdir(path):
            if fn.endswith('.pdb'):
                self.assertTrue(os.access(path+fn, os.F_OK))
                m = ModernaStructure('file', path+fn, 'A')
                print fn
                self.assertTrue(m['1'])
                self.assertTrue(m['2'])
                
            

if __name__ == '__main__':
    main()
