#!/usr/bin/env python

from unittest import TestCase, main, skip
from moderna.ModernaStructure import ModernaStructure
from moderna.sequence.ModernaSequence import Sequence
from moderna.Constants import LIR_DIRECTORY_PATH
from moderna.tests.test_data import *
import os

# from moderna.pycogent.rnaview_app import RnaView
# from moderna.pycogent.rnaview import RnaviewParser


class RnaViewTests(TestCase):

    @skip('requires rnaview installation')
    def test_calc_bp(self):
        """Should calculate all base pairs."""
        for fn in os.listdir(LIR_DIRECTORY_PATH):
            if not fn.endswith('.pdb'):continue
            pdb_file = LIR_DIRECTORY_PATH+fn
            print(pdb_file)
            #copy to temporary file
            tmp_file = '/tmp/rna.pdb'
            open(tmp_file, 'w').write(open(pdb_file).read())
            rna_view = RnaView()
            rna_view_result = rna_view(tmp_file)
            # parse result
            bpairs = rna_view_result['base_pairs']
            bp_dict = RnaviewParser(bpairs, strict=False)
            for bp in bp_dict['BP']:
                pass

    @skip('no assert found')
    def test_calc_bp_on_models(self):
        """Should calculate all base pairs."""
        for dname in os.listdir(MODEL_PATH):
            if not os.path.isdir(MODEL_PATH+dname): continue
            for tname in os.listdir(MODEL_PATH+dname):
                path = MODEL_PATH + dname + os.sep + tname + os.sep
                if not os.path.isdir(path): continue
                for fn in os.listdir(path):
                    if not fn.endswith('.pdb'):continue
                    pdb_file = path+fn
                    print(pdb_file)
                    #copy to temporary file
                    tmp_file = '/tmp/rna.pdb'
                    open(tmp_file, 'w').write(open(pdb_file).read())
                    rna_view = RnaView()
                    rna_view_result = rna_view(tmp_file)
                    # parse result
                    bpairs = rna_view_result['base_pairs']
                    bp_dict = RnaviewParser(bpairs, strict=False)
                    bplist = [bp for bp in bp_dict['BP']]

if __name__ == '__main__':
    main()
    
