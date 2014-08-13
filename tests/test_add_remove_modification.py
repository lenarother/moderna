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
from moderna.ModernaResidue import ModernaResidue
from moderna.analyze.BaseRecognizer import BaseRecognizer
from moderna.ModernaAlphabet import Alphabet
from moderna.Errors import ModernaResidueError
from moderna import load_model
import re
from test_data import *

class RemoveModificationTests(TestCase):
    """Makes sure modifications can be removed."""

    def test_remove(self):
        struc = ModernaStructure('file',MINI_TEMPLATE)
        r10 = struc['10']
        self.assertEqual(BaseRecognizer().identify_resi(r10),'m2G')
        r10.remove_modification()
        self.assertEqual(BaseRecognizer().identify_resi(r10),'G')

    def test_remove_deoxy(self):
        struc = ModernaStructure('file',DNA_WITH_MISMATCH,'E')
        r10 = struc['10']
        self.assertEqual(BaseRecognizer().identify_resi(r10),'dG')
        r10.remove_modification()
        self.assertEqual(BaseRecognizer().identify_resi(r10),'G')

    def test_remove_modification_empty(self):
        """Raises an Exception if base is not modified."""
        struc = ModernaStructure('file',MINI_TEMPLATE)
        self.assertRaises(ModernaResidueError,struc['11'].remove_modification)
        
    
class AddModificationTests(TestCase):
    """
    Makes sure modifications can be added to
    the four standard bases A,G,C,U in ModernaResidues.
    """
    def setUp(self):
        """Loads the A residue to start with."""
        self.struc = ModernaStructure('file',A_RESIDUE)
        self.adenosine = self.struc['1']

    def test_add_to_a(self):
        """Add modification to A."""
        self.adenosine.add_modification('m1A')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'m1A')

    def test_add_to_g(self):
        """Add modification to G."""
        self.adenosine.exchange_base('G')
        self.adenosine.add_modification('m1G')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'m1G')

    def test_add_to_u(self):
        """Add modification to U."""
        self.adenosine.exchange_base('U')
        self.adenosine.add_modification('Y')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'Y')

    def test_add_to_wrong_base(self):
        """Add modification to A that belongs to G should not work."""
        self.adenosine.add_modification('m1G')
        self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'m1G')
        atoms_m1G=["C1'",'C2',"C2'","C3'",'C4',"C4'",'C5',"C5'",'C6','C8','CM1','N1','N2','N3','N7','N9', "O2'","O3'","O4'","O5'",'O6','OP1','OP2','P']
        atoms=[at.name.strip() for at in self.adenosine.child_list]
        atoms.sort()
        self.assertEqual(atoms_m1G,atoms)
        
    def test_add_to_unk(self):
        """Should be possible to add modification to unknown residue when ribose is present"""
        m=load_model(PDB_UNK)
        for resi in m:
            resi.add_modification('m1G')
            self.assertEqual(resi.long_abbrev, 'm1G')
            
        
    def test_all(self):
        """Adding should work for all modifications."""
        a = Alphabet()
        br = BaseRecognizer()
        not_working = []
        errors = []
        EXCLUDED = ['A','G','C','U', 
                    '?A','?G','?C','?U',# exclude unknown
                    'X','?X','Xm', 'x',
                    'preQ0base','Qbase','preQ1base',
                    'galQtRNA',  # indistinguishable from ManQtRNA
                    '-', '_',
                    'yW-58','yW-72','yW-86','m8A','fa7d7G', # new in Modomics 2009, not yet in ModeRNA.
                    'm7Gpp_cap',  # not implemented yet
                    ]
        SYNONYMS = {'m42C':'m44C','m42Cm':'m44Cm','m62A':'m66A','m62Am':'m66Am'}
        for k in a:
            if k not in EXCLUDED and a[k].category not in ['unknown', 'standard', 'ligand', 'synthetic', 'stereoisomer', 'insertion', 'missing', ' ']:
                struc = ModernaStructure('file',A_RESIDUE)
                r = struc['1']
                try:
                    r.add_modification(k)
                    #print k, 'added'
                    #struc.write_pdb_file('dummies/'+k+'.pdb')
                    right = SYNONYMS.get(k,k)
                    if br.identify_resi(r) != right: 
                        not_working.append(k+','+br.identify_resi(r))
                        # write file for checking
                        struc.write_pdb_file('dummies/'+k+'.pdb')
                except ModernaResidueError:
                    raise
                    errors.append(k)
                                                       
        if not_working or errors:
            print '\nTest failed for modifications.'
            print 'Different base was recognized:'
            print ', '.join(not_working)
            print 'ERROR occured:'
            print ', '.join(errors)
         
        self.assertEqual(len(not_working)+len(errors),0)


class ExchangeModificationTests(TestCase):

    def setUp(self):
        """Loads the A residue to start with."""
        self.struc = ModernaStructure('file',A_RESIDUE)
        self.adenosine = self.struc['1']

    def test_mods_sanity(self):
        """Adding and removing many times should work as well."""
        #for mod in ['m1A','m66A','Am','t6A']:
        for mod in ['m1A','m6Am','Am','t6A']:
            # there is no modification named m66A. There is m6Am 
            self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),'A')
            self.adenosine.add_modification(mod)
            self.assertEqual(BaseRecognizer().identify_resi(self.adenosine),mod)
            self.adenosine.remove_modification()

    def test_dna_exchange(self):
        """All combinations of DNA->DNA exchanges should work."""
        bases = ['dT','dA','dG','dC']
        br = BaseRecognizer()
        r = self.adenosine
        for b1 in bases:
            r.add_modification(b1)
            self.assertEqual(br.identify_resi(r),b1)
            for b2 in bases:
                r.remove_modification()
                r.add_modification(b2)
                self.assertEqual(br.identify_resi(r),b2)

    

if __name__ == '__main__':
    main()
  
