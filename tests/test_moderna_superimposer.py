#!/usr/bin/env python
#
# test_moderna_superimposer.py
#
# unit tests for ModernaStructure.ModernaSuperimposer
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
from Bio.PDB import PDBParser
from moderna.ModernaSuperimposer import ModernaSuperimposer
from test_data import *
from moderna.Constants import BASE_PATH

A_SUPERIMPOSED_ON_C = {
'N1' :[ -3.89434432983, 8.9402217865, 46.3533706665],
'C2' :[ -4.63287353516, 9.30731010437, 45.3152999878],
'N3' :[ -4.31881713867, 9.30059432983, 44.0326004028],
'C4' :[ -3.04969406128, 8.84617042542, 43.8791046143],
'C5' :[ -2.18257522583, 8.44131088257, 44.8567810059],
'C6' :[ -2.64583969116, 8.48957633972, 46.1708526611],
'N6' :[ -1.920753479  , 8.07509803772, 47.2034835815],
'N7' :[ -0.969402313232, 8.05079460144, 44.3101501465],
'C8' :[ -1.13603591919, 8.22616958618, 43.0326538086],
'N9' :[ -2.37227630615, 8.70265960693, 42.6942443848]
}



class SuperimposerTests(TestCase):
    """
    Tests the functionality of ModernaSuperimposer as in the class
    description.
    """
    def test_superimposition(self):
        """
        """
        a = PDBParser().get_structure('A',BASE_PATH+'A.ent')[0]['C'][(' ',54,' ')]
        c = PDBParser().get_structure('C',BASE_PATH+'C.ent')[0]['C'][(' ',54,' ')] 
        original_coord = a['N6'].coord[1]
        ModernaSuperimposer([c['N1'],c['C2'],c['C6']],[a['N9'],a['C4'],a['C8']],a.child_list)
       
        self.failIfEqual(a['N6'].coord[1],original_coord)


    def test_superimposition2(self):
        """
        """
        a = PDBParser().get_structure('A',BASE_PATH+'A.ent')[0]['C'][(' ',54,' ')]
        c = PDBParser().get_structure('C',BASE_PATH+'C.ent')[0]['C'][(' ',54,' ')] 
        #original_coord = a['N6'].coord[1]
        ModernaSuperimposer([c['N1'],c['C2'],c['C6']],[a['N9'],a['C4'],a['C8']],a.child_list)
        for at in a:
            self.assertAlmostEqual(at.coord[0], A_SUPERIMPOSED_ON_C[at.id.strip()][0],5)
            self.assertAlmostEqual(at.coord[1], A_SUPERIMPOSED_ON_C[at.id.strip()][1],5)
            self.assertAlmostEqual(at.coord[2], A_SUPERIMPOSED_ON_C[at.id.strip()][2],5)

    def test_rmsd(self):
        """The RMS of the superposition should be close to manual superimp. in PyMOL."""
        a = PDBParser().get_structure('A',BASE_PATH+'A.ent')[0]['C'][(' ',54,' ')]
        g = PDBParser().get_structure('G',BASE_PATH+'G.ent')[0]['C'][(' ',54,' ')]
        m = ModernaSuperimposer([a['C4'],a['N9'],a['C8']],[g['C4'],g['N9'],g['C8']],a.child_list)
        self.assertTrue(m.rmsd < 0.004)



if __name__ == '__main__':
    main()
    
        
