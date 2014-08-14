#!/usr/bin/env python
#
# ModelingTask.py
#
# Operations that build the model.
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


class ModelingTask(object):
    """Abstract superclass"""
    def build(self, template, model):
        pass
        
class SingleResidueTask(ModelingTask):
    def __init__(self, residue, new_number=None):
        self.residue = ModeRNAResidue(residue)
        self.number = new_number
        
class CopyResidueTask(SingleResidueTask):
    
    def build(self, template, model):
        """Copies a single residue to the model."""
        B_FACTOR_COPY = True
        #TODO: add Bfactor option
        num = self.number or self.residue.identifier
        self.residue.set_bfactor(B_FACTOR_COPY)
        model.add_residue(self.residue, num, strict = strict)       
        log.write_message('Residue %s: residue copied from template to model.' %num)
        
        
