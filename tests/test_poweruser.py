#!/usr/bin/env python
#
# test_poweruser.py
#
# unit tests for base exchange functionality
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

from unittest import main, TestCase
from moderna.ModernaStructure import ModernaResidue,ModernaStructure
import os

RNADB_PATH = '../../data/rnaDB05/'

MISC_PATH = '../../../suffix/sourceforge/trunk/data/'

class PoweruserTests(TestCase):
    """
    Checks how Moderna handles a wide range of structures.
    """
    def test_load_rnadb(self):
        """Should load all of Richardsons RNADB2005 structures."""
        self.check_struc_dir(RNADB_PATH)

    def test_load_hires(self):
        """Should load Tomek Osinskis Hires dataset."""
        self.check_struc_dir(MISC_PATH+'hi-res_cleaned/')

    def test_load_trna(self):
        """Should load all of Rother/Bauers tRNA structures."""
        self.check_struc_dir(MISC_PATH+'trna\\')
        self.check_struc_dir(MISC_PATH+'trna_curated\\')

    def test_load_risc(self):
        self.check_struc_dir(MISC_PATH+'risc\\')

    def test_load_(self):
        self.check_struc_dir()

    def check_struc_dir(self,path):
        loaded = 0
        succeeded = 0
        for fn in os.listdir(path):
            if fn[-4:].lower() not in ['.ent','.pdb']: continue
            loaded += 1
            try:
                s = ModernaStructure('file',path+fn)
                succeeded += 1
            except Exception,e:
                print 'failed to load:',path+fn
                print e
        self.assertEqual(succeeded,loaded)


if __name__ == '__main__':
    main()
  
