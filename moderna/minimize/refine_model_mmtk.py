#!/usr/bin/env python
#
# refine_model_mmtk.py
#
# optimizes models with MMTK

import os, sys, tempfile
from optparse import OptionParser

from moderna.moderna import *
from moderna.RNAModel import RnaModel
from moderna.ModernaStructure import ModernaStructure
from moderna.ModernaSuperimposer import ModernaSuperimposer
from moderna.Constants import BACKBONE_ATOMS
#from moderna.dihedral import angle
from Bio.PDB.Vector import Vector, calc_angle, calc_dihedral

from MMTK import *
from MMTK.ForceFields import Amber94ForceField
from MMTK.Minimization import SteepestDescentMinimizer
from MMTK.Minimization import ConjugateGradientMinimizer
from MMTK.ForceFields.Restraints import HarmonicDistanceRestraint, HarmonicAngleRestraint, HarmonicDihedralRestraint
from MMTK.PDB import PDBConfiguration
from MMTK.NucleicAcids import NucleotideChain
from MMTK.Visualization import view

TEMPFILE_BASE = 'temp_mmtk_'

class OptimizationError(Exception): pass
class NameReplacementError(OptimizationError): pass

class NameReplacement:
    """Class for replacing atom and residue names."""
    atoms_amber = [' O1P', ' O2P', ' O3P', ' C1*', ' C2*', ' C3*', ' C4*', ' C5*', ' O2*', ' O3*', ' O4*', ' O5*']
    atoms_moderna = [" OP1", " OP2", " OP3", " C1'", " C2'", " C3'", " C4'", " C5'", " O2'", " O3'", " O4'", " O5'"]
    residues_amber = [' RA', ' RC', ' RG', ' RU']
    residues_moderna = ['  A', '  C', '  G', '  U']

    def __init__(self, mode):
        if mode == 'back':
            self.atom_dict = dict(zip(self.atoms_amber, self.atoms_moderna))
            self.resi_dict = dict(zip(self.residues_amber, self.residues_moderna))
        elif mode == 'forth':
            self.atom_dict = dict(zip(self.atoms_moderna, self.atoms_amber))
            self.resi_dict = dict(zip(self.residues_moderna, self.residues_amber))
        else:
            raise NameReplacementError("Invalid options for name replacement")

    def replace_residue_name(self, name):
        return self.resi_dict.get(name, name)

    def replace_atom_name(self, name):
        return self.atom_dict.get(name, name)


BACK_SET = NameReplacement('back')
FORTH_SET = NameReplacement('forth')

class ModelMinimization:
    """
    Main class for minimizing models
    """
    def __init__(self, options, temp_path):
        self.options = options
        # structure
        print 'Loading model from %s, chain %s'%(options.input_file, options.chain_name)
        self.model = load_model(options.input_file, options.chain_name)
        self.model_passive = None
        self.stem5,  self.stem3 = None, None
        self.sequence_before = self.model.get_sequence()
        self.model_resnumbers = []
        # parameters
        self.residues = None
        if options.residues:
            self.residues = options.residues.split('-')
        self.cycles = int(options.cycles)
        self.output_name = options.output_file
        self.modifications = []
        self.restraints = []
        self.mmtk_output_file = 'mmtk_output.pdb'
        self.temp_pdb_file = temp_path

    def model_to_tempfile(self):
        """writes model to temporary file."""
        self.model.write_pdb_file(self.temp_pdb_file)
        self.model = None

    def tempfile_to_model(self):
        """Reads model from temporary file."""
        self.model = load_model(self.temp_pdb_file, self.options.chain_name)

    def extract_region(self):
        self.model.renumber_chain()
        """Cuts out the residues for optimization and keeps the rest in memory."""
        if self.residues:
            print '\nExtracting the region to be optimized (%s-%s) from the model.'%(self.residues[0], self.residues[1])
            print '\t.. total residues in model                    : %4i'%len(self.model)
            self.model_passive = self.model.find_residues_not_in_range(self.residues[0], self.residues[1])
            print '\t.. residues not participating in optimization : %4i (+2 for superposition afterwards)'%(len(self.model_passive)-2)
            print '\t   '+','.join([r.identifier for r in self.model_passive])
            self.model = RnaModel(None, None, self.options.chain_name, 'residues', self.model[self.residues[0]:self.residues[1]])
            # memorize residue numbers because MMTK screws them
            self.model_resnumbers = [r.identifier for r in self.model]
            print '\t.. residues participating in optimization     : %4i'%len(self.model)
            print '\t   '+','.join([r.identifier for r in self.model])
            # keeps duplicate versions of the residues indicated by self.residues,
            # so they can be superimposed later
            self.stem5,  self.stem3 = self.model[self.residues[0]], self.model[self.residues[1]]
        else:
            self.struc_to_optimize = self.model

    def merge_region(self):
        """Merges the optimized residues with the memorized ones."""
        residues = [r for r in self.model]
        if self.model_passive:
            # apply old residue numbers
            print '\nRestoring original numeration in the optimized region.'
            for num, resi in zip(self.model_resnumbers, residues):
                #print resi, num
                #self.model.renumber_residue(resi.identifier, num)
                resi.change_number(num)
            self.model = ModernaStructure('residues',residues)
            # do the superposition
            print '\nSuperimposing the optimized part onto the rest of the model.'
            all_atoms = self.model.get_all_atoms()
            sup = ModernaSuperimposer(moved_atoms = all_atoms)
            sup.get_atoms([self.stem5, self.stem3], BACKBONE_ATOMS, 'fixed')
            resi5, resi3 = self.model[self.residues[0]], self.model[self.residues[1]]
            sup.get_atoms([resi5, resi3], BACKBONE_ATOMS, 'moved')
            sup.superimpose()

            # merge residues
            print '\nMerging the optimized part with the rest of the model.'
            resi = [r for r in self.model]
            for r in self.model_passive:
                if r.identifier not in self.residues:
                    resi.append(r)
            self.model = RnaModel(None, None, self.options.chain_name, 'residues', resi)


    def remove_modifications(self):
        """Memorizes modified nucleotides for later, and then removes them."""
        print '\nRemoving modified bases from the model.'
        for resi in self.model.get_modified_residues():
            resi=self.model[resi]
            modi = (resi.identifier, resi.long_abbrev)
            print '\t.. %s at position %s'%modi
            self.modifications.append(modi)
        self.model.remove_all_modifications()
        # check for O3 atoms
        for resi in self.model:
            if resi.child_dict.has_key(' O3P'):
                raise OptimizationError('O3P atoms cannot be handled by MMTK!')


    def change_atom_names(self, forth=False, back=False):
        """Replaces atom and residue names so that MMTK can iterpret them."""
        print '\nReplacing names in temporary PDB file (%s)'%((forth and 'ModeRNA->MMTK') or (back and 'MMTK->ModeRNA'))
        if forth: nameset = FORTH_SET
        elif back: nameset = BACK_SET
        else: raise OptimizationError("Invalid options for changing atom names.")
        # replace residue and atom names
        out = []
        for line in open(self.temp_pdb_file):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resiname = nameset.replace_residue_name(line[17:20])
                atomname = nameset.replace_atom_name(line[12:16])
                line = line[:12]+atomname+line[16:17]+resiname+line[20:]
            if len(line)>13 and line[13]!='H': # discard hydrogens
                out.append(line)
        open(self.temp_pdb_file, 'w').writelines(out)


    def create_restraints(self, chain):
        """If a residue range was specified, these residues are cut out, and """
        if self.residues:
            print '\t.. building restraints'
            first_resi = chain[0]
            last_resi = chain[-1]

            # distances between the two terminal residues
            dist_po = (self.stem5['P'] - self.stem3["O3'"])/10.0
            restraints = HarmonicDistanceRestraint(first_resi.phosphate.P, last_resi.sugar.O_3, dist_po, 10000.)
            dist_co = (self.stem5["O5'"] - self.stem3["O3'"])/10.0
            restraints += HarmonicDistanceRestraint(first_resi.sugar.O_5, last_resi.sugar.C_3, dist_co, 10000.)
            dist_cc = (self.stem5["C5'"] - self.stem3["C4'"])/10.0
            restraints += HarmonicDistanceRestraint(first_resi.sugar.C_5, last_resi.sugar.C_4, dist_cc, 10000.)

            # angles between the two terminal residues
            angle_cpo = calc_angle(self.stem5["C4'"].get_vector(), self.stem3["P"].get_vector(), self.stem3["O3'"].get_vector())
            #angle_cpo = calc_angle(self.stem3["P"].coord-self.stem5["C4'"].coord, self.stem3["P"].coord-self.stem3["O3'"].coord)
            restraints += HarmonicAngleRestraint(first_resi.sugar.C_4, last_resi.phosphate.P, last_resi.sugar.O_3, angle_cpo, 100000.)
            angle_pco1 = calc_angle(self.stem5["P"].get_vector(), self.stem3["C4'"].get_vector(), self.stem3["O3'"].get_vector())
            restraints += HarmonicAngleRestraint(first_resi.phosphate.P, last_resi.sugar.C_4, last_resi.sugar.O_3, angle_pco1, 10000.)
            angle_pco2 = calc_angle(self.stem5["P"].get_vector(), self.stem5["C4'"].get_vector(), self.stem3["O3'"].get_vector())
            restraints += HarmonicAngleRestraint(first_resi.phosphate.P, first_resi.sugar.C_4, last_resi.sugar.O_3, angle_pco2, 10000.)
            angle_pcc = calc_angle(self.stem5["P"].get_vector(), self.stem5["C4'"].get_vector(), self.stem3["C4'"].get_vector())
            restraints += HarmonicAngleRestraint(first_resi.phosphate.P, first_resi.sugar.C_4, last_resi.sugar.C_4, angle_pcc, 10000.)
            # torsion between the two terminal residues
            # THESE PRODUCE SEGMENTATION FAULTS - THE OTHERS MUST SUFFICE
            # restraints += HarmonicDihedralRestraint(first_resi.phosphate.P, first_resi.sugar.C_4, last_resi.phosphate.P, last_resi.sugar.C_4, 1.0, 10.)
            #restraints += HarmonicDihedralRestraint(first_resi.phosphate.P, first_resi.sugar.C_4, last_resi.sugar.C_4, last_resi.sugar.O_3, 0.0, 100000.)
            #restraints += HarmonicDihedralRestraint(first_resi.sugar.C_4, first_resi.sugar.O_3, last_resi.phosphate.P, last_resi.sugar.C_4, 0.0, 100000.)
            #restraints += HarmonicDihedralRestraint(first_resi.sugar.C_4, first_resi.sugar.O_3, last_resi.sugar.C_4, last_resi.sugar.O_3, 0.0, 100000.)

            return restraints

    def optimize(self):
        """Run the optimization with MMTK"""
        print '\n-------------------------------------------------------------------------------'
        print '\n MMTK Optimization starts...'

        print '\t.. building universe'
        configuration = PDBConfiguration(self.temp_pdb_file)
        # Construct the nucleotide chain object. This also constructs positions
        # for the missing hydrogens, using geometrical criteria.
        chain = configuration.createNucleotideChains()[0]
        universe = InfiniteUniverse()
        universe.addObject(chain)

        restraints = self.create_restraints(chain)

        # define force field
        print '\t.. setting up force field'
        if restraints:
            ff = Amber94ForceField() + restraints
        else:
            ff = Amber94ForceField()
        universe.setForceField(ff)

        # do the minimization
        print '\t.. starting minimization with %i cycles'%self.cycles
        minimizer = ConjugateGradientMinimizer(universe)
        minimizer(steps = self.cycles)
        # write the intermediate output
        print '\t.. writing MMTK output to %s'% self.temp_pdb_file
        if self.model_passive:
            print '\t   (please note that MMTK applies a different numeration of residues.\n\t    The original one will be restored in the final output).'
        universe.writeToFile(self.mmtk_output_file)
        open(self.temp_pdb_file, 'w').write(open(self.mmtk_output_file).read())
        print '\n-------------------------------------------------------------------------------'

    def add_modifications(self):
        """Restores modifications on the model."""
        print "\nAdding modifications back to the model:"
        for position, modif in self.modifications:
            print "\t.. adding %s in position %s"%(modif, position)
            self.model[position].add_modification(modif)
        # check if sequence stayed the same
        print '\nSequence before optimization'
        print self.sequence_before
        print 'Sequence after optimization'
        print self.model.get_sequence()

    def write_result(self):
        """Creates the output file."""
        print "\nWriting final output to: %s"%self.output_name
        self.model.write_pdb_file(self.output_name)



def parse_options():
    """
-------------------------------------------------------------------------------
    ModeRNA (C) 2009 by Magdalena Musielak, Kristian Rother et al.

    Minimizer for RNA models using MMTK.

    support: krother@genesilico.pl, t.puton@amu.edu.pl
-------------------------------------------------------------------------------
    """
    usage = """usage: %prog [options] arg

Examples:
refine_model_mmtk.py -m starting_structure.pdb -c A -y 250 -o optimized_structure.pdb
refine_model_mmtk.py -m starting_structure.pdb -c A -y 250 -r 20-50 -o optimized_region.pdb"""
    print parse_options.__doc__
    parser = OptionParser(usage, version="%prog 1.0")

    parser.add_option("-m", "--model", dest="input_file", help="PDB file with the model to optimize")
    parser.add_option("-c", "--chain", dest="chain_name", default="A", help="Chain to read (default is A)")
    parser.add_option("-o", "--output", dest="output_file", help="PDB file to which output is written")
    parser.add_option("-y", "--cycles", dest="cycles", help="number of optimization cycles")
    parser.add_option("-r", "--residues", dest="residues", help="range of the model to optimize (e.g. '10-20').")
    parser.add_option("-t", "--remove_temp", dest="remove_temp",  action="store_true", default=False, help="remove temporary PDB file")

    (options, args) = parser.parse_args()
    return options

def validate_options(options):

    if not (options.input_file and options.output_file and options.cycles):
        raise OptimizationError('Invalid input. Type -h to see all available options.\n')
    if not os.access(options.input_file, os.F_OK):
        raise OptimizationError("Input file '%s' does not exist.\n"%options.input_file)
    try:
        io = int(options.cycles)
    except:
        raise OptimizationError("Number of cycles must be a number (1-10000).\n")
    if not 1<io<10000:
        raise OptimizationError("Number of cycles must be a number (1-10000).\n")
    if options.residues:
        sp = options.residues.split('-')
        if len(sp)!=2:
            raise OptimizationError("Residue numbers must specify a range, e.g. '10-20'.\n")


def main(options):
    """Performs the optimization."""
    temp_fd, temp_path = \
        tempfile.mkstemp(prefix=TEMPFILE_BASE, suffix=".pdb", dir=os.getcwd())
    mm = ModelMinimization(options, temp_path)
    mm.extract_region()
    mm.remove_modifications()
    mm.model_to_tempfile()
    # --- now the program is working on a PDB file instead of a ModeRNA object
    mm.change_atom_names(forth=True)
    mm.optimize()
    mm.change_atom_names(back=True)
    # --- end part working on file
    mm.tempfile_to_model()
    mm.merge_region()
    mm.add_modifications()
    mm.write_result()
    if options.remove_temp:
        print 'Removing temporary PDB file: %s' % temp_path
        os.close(temp_fd)
        os.unlink(temp_path)
    else:
        print 'Temporary PDB file has not been removed: %s' % temp_path

if __name__ == "__main__":
    try:
        options=parse_options()
        validate_options(options)
        main(options)
    except OptimizationError, e:
        sys.stderr.write(str(e))

