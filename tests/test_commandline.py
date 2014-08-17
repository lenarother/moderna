#!/usr/bin/env python
#
# test_moderna.py
#
# unit tests for moderna interface
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
from moderna.Template import Template
from moderna.sequence.ModernaSequence import Sequence
from test_data import *
from moderna.Constants import MODULE_PATH
import os
import tempfile

this_dir = os.getcwd()

class  CommandlineTests(TestCase):
    """
    Command-line-interface tests
    """
    def setUp(self):
        self.tmp = tempfile.mktemp()
        self.tmp2 = tempfile.mktemp()
        
    def tearDown(self):
        if os.access(self.tmp,os.F_OK): os.remove(self.tmp)
        if os.access(self.tmp2,os.F_OK): os.remove(self.tmp2)
        
    def test_build_model(self):
        cmd = "python %smoderna.py -t %s -a %s -o %s > %s"%(MODULE_PATH+os.sep, MINI_TEMPLATE,MINI_ALIGNMENT_FILE,self.tmp, self.tmp2)
        os.system(cmd)
        t = Template(self.tmp)
        self.assertEqual(t.get_sequence(),Sequence('ACUGUGAYUA[UACCU#PG'))
        
    def test_insert_modification(self):
        cmd = "python %smoderna.py -s %s -o %s -m m22G -p 3 > %s"%(MODULE_PATH+os.sep, MINI_TEMPLATE, self.tmp, self.tmp2)
        os.system(cmd)
        t = Template(self.tmp)
        self.assertEqual(t.get_sequence(),Sequence('GCRGAUUUALCUCAG'))

    def test_insert_modi_chain(self):
        cmd = "python %smoderna.py -s %s -c D -o %s -m m22G -p 38 > %s"%(MODULE_PATH+os.sep, RNA_HAIRPIN, self.tmp, self.tmp2)
        os.system(cmd)
        t = Template(self.tmp, template_chain_name='D')
        self.assertEqual(t.get_sequence(),Sequence('CUGCCUQUC/RGC'))

    def test_get_sequence(self):
        cmd = "python %smoderna.py -t %s > %s"%(MODULE_PATH+os.sep, MINI_TEMPLATE, self.tmp)
        os.system(cmd)
        seq = open(self.tmp).read().strip()
        self.assertEqual(seq, "TEMPLATE SEQUENCE:\nGCGGAUUUALCUCAG")
        
    def test_examine_template(self):
        cmd = "python %smoderna.py -t %s -e > %s"%(MODULE_PATH+os.sep, NASTY_PDB, self.tmp)
        os.system(cmd)
        report = open(self.tmp).read().strip()
        self.assertTrue(report.find("Water residues present in the chain")>-1)
        
    def test_examine_template(self):
        cmd = "python %smoderna.py -t %s -l > %s"%(MODULE_PATH+os.sep, NASTY_PDB, self.tmp)
        os.system(cmd)
        report = open(self.tmp).read().strip()
        self.assertFalse(report.find("._._._")>-1)
        


if __name__ == '__main__':
    main()
    
