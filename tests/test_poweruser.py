#!/usr/bin/env python


from unittest import main, TestCase, skip
from moderna.RNAResidue import RNAResidue
from moderna.ModernaStructure import ModernaStructure
import os

RNADB_PATH = '../../data/rnaDB05/'

MISC_PATH = '../../../suffix/sourceforge/trunk/data/'


class PoweruserTests(TestCase):
    """
    Checks how Moderna handles a wide range of structures.
    """
    @skip('requires customizing of test runner or tag')
    def test_load_rnadb(self):
        """Should load all of Richardsons RNADB2005 structures."""
        self.check_struc_dir(RNADB_PATH)

    @skip('requires customizing of test runner or tag')
    def test_load_hires(self):
        """Should load Tomek Osinskis Hires dataset."""
        self.check_struc_dir(MISC_PATH+'hi-res_cleaned/')

    @skip('requires customizing of test runner or tag')
    def test_load_trna(self):
        """Should load all of Rother/Bauers tRNA structures."""
        self.check_struc_dir(MISC_PATH+'trna\\')
        self.check_struc_dir(MISC_PATH+'trna_curated\\')

    @skip('requires customizing of test runner or tag')
    def test_load_risc(self):
        self.check_struc_dir(MISC_PATH+'risc\\')

    @skip('requires customizing of test runner or tag')
    def test_load_(self):
        self.check_struc_dir()

    def check_struc_dir(self,path):
        loaded = 0
        succeeded = 0
        for fn in os.listdir(path):
            if fn[-4:].lower() not in ['.ent','.pdb']:
                continue
            loaded += 1
            try:
                s = ModernaStructure('file', path + fn)
                succeeded += 1
            except Exception:
                print('failed to load:', path + fn)
        self.assertEqual(succeeded, loaded)


if __name__ == '__main__':
    main()
  
