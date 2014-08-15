#!/usr/bin/env python
#
# test_alphabet.py
#
# creates the G,C,U standard bases from a sample A.
#
__author__ = "Magdalena Musielak, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Magdalena Musielak"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Prototype"

from moderna.ModernaStructure import ModernaResidue,ModernaStructure
from moderna.analyze.BaseRecognizer import BaseRecognizer

NUC_PATH = 'standard_nucleotides/'
A_RESIDUE = NUC_PATH+'a.ent'

if __name__ == '__main__':
    # Loads the A residue to start with.
    struc = ModernaStructure('file',A_RESIDUE)
    ade = struc[1]

    # exchange to G
    ade.exchange_base('G')
    assert(BaseRecognizer().identify_resi(ade),'G')
    struc.write_pdb_file(NUC_PATH+'g.ent')
    # exchange to C
    ade.exchange_base('C')
    assert(BaseRecognizer().identify_resi(ade),'C')
    struc.write_pdb_file(NUC_PATH+'c.ent')
    # exchange to U
    ade.exchange_base('U')
    assert(BaseRecognizer().identify_resi(ade),'U')
    struc.write_pdb_file(NUC_PATH+'u.ent')

