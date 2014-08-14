#!/usr/bin/env python
#
# setup.py
#
# Generate Windows executable version of ModeRNA.
# 
# http://iimcb.genesilico.pl/moderna/ 
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "1.7.1"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"


from distutils.core import setup
import os
import os.path

#
# setup script to generate Windows executable only!
#
if os.sep!='/':
    import py2exe

def all_from_dir(path):
    if path[-1] != os.sep:
        path += os.sep
    files = [path+f for f in os.listdir(path) if not f=='.svn']
    return path, files

file_list = [
              ('.',['LICENSE_GPL.TXT','README.TXT','RELEASE_NOTES.TXT']),
              ('data'+os.sep,[
                  'data'+os.sep+'LIR_fragments.lib',
                  'data'+os.sep+'rnaDB05_list.txt',
                  'data'+os.sep+'modification_names_table',
                  'data'+os.sep+'modification_topologies.txt',
                  'data'+os.sep+'modifications',
                  'data'+os.sep+'helix.pdb',
                  'data'+os.sep+'suite_clusters.txt',
                  'data'+os.sep+'pair_fragment.pdb',
                  'data'+os.sep+'single_strand.pdb',
                  'data'+os.sep+'phosphate_group.pdb',
                  'data'+os.sep+'IsostericityMatrices.txt',
              ]),
               all_from_dir('data'+os.sep+'modification_fragments'),
               all_from_dir('data'+os.sep+'rnaDB05'),
               all_from_dir('data'+os.sep+'standard_bases'),
               all_from_dir('data'+os.sep+'base_pairs'),
    ]

setup(
        name = "ModeRNA",
        console=['moderna.py'],
        version = "1.6.0",
        author="Magdalena Musielak, Kristian Rother, Tomasz Puton, Janusz M. Bujnicki",
        author_email="mmusiel@genesilico.pl",
        url="http://iimcb.genesilico.pl/moderna", 
        py_modules  = [
         'AlignmentMatcher', 
         'BackboneBuilder','BaseRecognizer','BasepairCalculator',
         'CheckPdb', 
         'ClashRecognizer','Constants','CoordBuilder',
         'decorators', 'Errors',
         'FCCDLoopCloser',
         'GeometryParameters','GeometryAnalyzer','GeometryStatistics',
         'HBondCalculator',
         'IsostericityMatrices','Isostericity', #'IsostericModelling',
         'LIRdb','LogFile','LIR',
         #'model_to_adun','madun',
         'ModernaAlignment','ModernaResidue','ModernaSequence',
         'ModernaAlphabet','ModernaFragment','ModernaStructure',
         'ModernaSuperimposer','ModernaToAdun',
         'ProcessPDB','PhosphateBuilder','RNAModel', 'PuckerCalculator', 
         'Renumerator', 'RNAAlignment','RNAChain', 'RNAResidue', 'RNASuites', 
         'StructureLibrary','SearchLIR','StackingCalculator',
         'SecstrucFragment', 
         'Template','moderna',
         '__init__','setup',
            ], 
        data_files=file_list)
        




