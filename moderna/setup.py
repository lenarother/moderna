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
import inspect
import Errors
#
# setup script to generate Windows executable only!
#
if os.sep!='/':
    import py2exe

MPATH = os.path.dirname(inspect.getfile(Errors))
print MPATH
DATA = MPATH + os.sep + 'data' + os.sep

def all_from_dir(path):
    path = os.path.abspath(path)
    if path[-1] != os.sep:
        path += os.sep
    files = [path + f for f in os.listdir(path) if not f=='.svn']
    return path, files

file_list = [
              ('.', ['LICENSE_GPL.TXT','README.TXT','RELEASE_NOTES.TXT']),
              (DATA, [
                  DATA + 'LIR_fragments.lib',
                  DATA + 'rnaDB05_list.txt',
                  DATA + 'modification_names_table',
                  DATA + 'modification_topologies.txt',
                  DATA + 'modifications',
                  DATA + 'helix.pdb',
                  DATA + 'suite_clusters.txt',
                  DATA + 'pair_fragment.pdb',
                  DATA + 'single_strand.pdb',
                  DATA + 'phosphate_group.pdb',
                  DATA + 'IsostericityMatrices.txt',
              ]),
               all_from_dir( 'modification_fragments'), # DATA +
               all_from_dir(DATA + 'rnaDB05'),
               all_from_dir(DATA + 'standard_bases'),
               all_from_dir(DATA + 'base_pairs'),
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
        




