
"""
Script for creating single-residue PDB structures for
all modifications in MODOMICS.
"""

from moderna.ModernaStructure import ModernaResidue,ModernaStructure
from modomics.modifications import get_modifications

NORMAL_BASE = 'standard_nucleotides/a.ent'
OUTPUT_PATH = 'modified_nucleotides/'

def create_modification_files():
    """
    Builds a model for each modification.
    """
    f=open(OUTPUT_PATH+'create_modifications_log','w')
    for mod_name in get_modifications():
        f.write(mod_name+'\n')
        # modifications_to_add=get_all_modifications_names()
        s=ModernaStructure('file',NORMAL_BASE)
        print s
        try:
            print s['1']
            s['1'].add_modification(mod_name)
            f.write('created structure for %s'%mod_name)
            f.write('\n'+'-'*30+'\n')
        except:
            f.write('no rule for %s.'%(mod_name))
            f.write('\n'+'-'*30+'\n')
        s.write_pdb_file(OUTPUT_PATH + mod_name+'.pdb')


if __name__ == '__main__':
    create_modification_files()
            
