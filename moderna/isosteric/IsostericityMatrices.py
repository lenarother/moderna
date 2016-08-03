#!/usr/bin/env python
"""
Module for modeling isosteric base pairs.
"""


from moderna.Constants import DATA_PATH
from moderna.util.LogFile import log


class IsostericityMatrices:
    """
    Isostericity matrices implementation in Python
    """
    def __init__(self):
        self.source = open(DATA_PATH+"IsostericityMatrices.txt","r")
        self.matrices = self.import_matrices()
        self.source.close()
        
    def import_matrices(self):
        """
        Reads data from source txt files and prepare the dictionary
        """
        matrices = {}
        for line in self.source:
            if line.startswith("Type: "):
                master_key = line.split(" ")[1]
                type_key = line.split(" ")[2].strip()
                if master_key not in matrices.keys():
                    matrices[master_key] = {}
                matrices[master_key][type_key] = {}
            elif line.startswith("\n"):
                continue
            else:
                data = line.split(": ")
                if data[1] == '\n':
                    matrices[master_key][type_key][data[0]] = None
                else:
                    matrices[master_key][type_key][data[0]] = float(data[1].strip()) 
        return matrices 
        
    def check_isostericity(self, bp1, bp2, interact_type, max_value=1.0):
        """
        Returns True if basepair1 is isosteric to basepair2 when interaction type is interact_type
        """
        try:
            result = self.matrices[bp1][interact_type][bp2] 
            log.write_message(bp1+"->"+bp2+" ("+interact_type+")")
            return result <= max_value and result is not None
        except:
            log.write_message("No information in IsostericityMatrices about: "+bp1+"->"+bp2+" ("+interact_type+")")
            return False
        
        
    def show_isosteric_bp(self, bp1, interact_type, max_value=1.0):
        """
        Returns a tuple with all base pairs isosteric to bp1 when interaction type is interact_type
        """
        pairs = self.matrices[bp1][interact_type]
        result = [bp for bp in pairs if pairs[bp] is not None and pairs[bp] <= max_value]
        return tuple(result)
