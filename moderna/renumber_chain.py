#!/usr/bin/env python
#
# renumber_chain.py
#
# Changes numbers in RNAStructures.
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

from LogFile import log
from Renumerator import NumberRenumerator
from Errors import RenumeratorError

def renumber_chain(chain, start_number='1'):
    """
    Renumbers an entire RNA chain object with continuous integer numbers.
    Starts from the given number (by default 1).

    Arguments:
    * start number (by default 1)
    """
    try:
        num = int(start_number)
    except ValueError: 
        raise RenumeratorError('Cannot start numeration with %s, requires number' % str(start_number)) 
      
    chain.sort_residues()
    numbers = map(str, range(num, num + len(chain)))
    renumber_all_residues(chain, numbers)        
    log.write_message('Chain was renumbered. The first residue is now %s.' % str(start_number))
        
        
def renumber_all_residues(chain, numbers):
    '''
    Changes the numbers of all residues.
    '''
    if len(chain) != len(numbers):
        raise RenumeratorError('amount of numbers (%i) does not match number of residues (%i)' % (len(numbers), len(chain))) 
    residues = list(chain)
    chain.remove_all_residues()
    for resi, num in zip(residues, numbers):
        chain.add_residue(resi, num)
    log.write_message('Chain was renumbered. New numbers are now %s.' % str(numbers))
      
      
def insert_gap(chain, after_number, num_residues):
    """
    Changes the numeration so that Renumerator can insert additional residue numbers.

    Arguments:
    * residue number, after which the gap will be inserted (None=at the beginning)
    * number of residues in the gap
    """
    chain.sort_residues()
    numbers = [resi.identifier for resi in chain]
    insert = [None] * num_residues    
    if after_number == None:
        new_numbers = insert + numbers
        idx = -1
    else:
        after_number = str(after_number)
        idx = numbers.index(after_number)
        new_numbers = numbers[:idx+1] + insert + numbers[idx+1:]
    renumerator = NumberRenumerator(new_numbers)
    new_numbers = renumerator.get_identifiers()
    new_numbers = new_numbers[:idx+1] + new_numbers[idx+1+num_residues:]
    renumber_all_residues(chain, new_numbers)

    log.write_message("%s-residue-long space inserted after residue %s." %\
        (str(num_residues), str(after_number)))

