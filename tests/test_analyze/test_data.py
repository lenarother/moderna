#!/usr/bin/env python
"""
Some typical test examples
"""

import os

MODULE_PATH,module_name = os.path.split(__file__)
if not MODULE_PATH:
    MODULE_PATH = os.getcwd()

TEST_DATA_PATH = MODULE_PATH + '/../test_data/'
#
# auxiliary files
#
TEST_MODIFICATION_TABLE = TEST_DATA_PATH+'modification_sample_table.txt'

#
# Nucleotide and modification examples
#

# single adenine
A_RESIDUE = TEST_DATA_PATH+'nucleotides/a.ent'
C_RESIDUE = TEST_DATA_PATH+'nucleotides/c.ent'

# with hydrogens, chain B
MOD_WITH_HYDROGENS =TEST_DATA_PATH+'nucleotides/mod_with_hydrogens.pdb'

# 2'3' phosphate diester group
THREEPRIME_PHOSPHATE = TEST_DATA_PATH+'nucleotides/23p.pdb'

# broken residues
A_TOO_SHORT = TEST_DATA_PATH+'nucleotides/a_too_short_glycosidic.ent'
UNRECOGNIZABLE_RESIDUE = TEST_DATA_PATH+'nucleotides/bad_residue.pdb'

# simple base pair
AC_BASEPAIR = TEST_DATA_PATH+'nucleotides/ac_basepair.pdb' # chain 0
DISTANT_AC_BASEPAIR = TEST_DATA_PATH+'nucleotides/ac_basepair_a498_c494_rr0011_distorted.ent'

# residue with just two H atoms
H_ONLY = TEST_DATA_PATH+'nucleotides/h_only_residue.pdb'

#
# RNA structures
#
RNA_PATH = TEST_DATA_PATH+'rna_structures/'

# 14-residue single stranded RNA
MINI_TEMPLATE = TEST_DATA_PATH+'rna_structures/minirna_14res_ss.pdb'
MINI_TEMPLATE_CHAIN_NAME = 'A'

# hires tRNA structure with modifications, chain A
RNA_1EHZ = TEST_DATA_PATH+'rna_structures/1ehz.ent'

# RNA arn035H from Richardson set
RICHARDSON_STRUC = TEST_DATA_PATH+'rna_structures/1rna.ent'

# some other RNA
RNA_HAIRPIN = TEST_DATA_PATH + 'rna_structures/rna_hairpin_fragment_D.pdb'
RNA_1QF6 = TEST_DATA_PATH+'rna_structures/1qf6.ent'
RNA_2OW8 = TEST_DATA_PATH+'rna_structures/2ow8.ent'
RNA_1B23 = TEST_DATA_PATH+'rna_structures/1b23_R.ent'
RNA_1C0A = TEST_DATA_PATH+'rna_structures/1c0a_B.ent'
RNA_2BTE =TEST_DATA_PATH+'rna_structures/2bte_B.ent'
RNA_1EXD = TEST_DATA_PATH+'rna_structures/1exd_B.ent'
RNA_1EFW = TEST_DATA_PATH+'rna_structures/1efw_D.ent'
JMB_TEMPLATE = TEST_DATA_PATH+'rna_structures/template_jmb.pdb'
MESSED_ATOMS = TEST_DATA_PATH+'rna_structures/atoms_messed_up.pdb'
# RNA with hydrogens
RNA_HYDRO = TEST_DATA_PATH+'rna_structures/rna_with_hydro_B.ent'
# RNA with two big strand breaks
DOUBLEGAP = TEST_DATA_PATH+'rna_structures/rna_2breaks.ent'
# with atoms and no header
RNA_ATOMSONLY = RNA_PATH + 'rna_atomsonly.ent'

CAP_EXAMPLE = RNA_PATH + 'mrna_cap5.pdb'

# PDB file with negative numbers
NEG_NUM_PDB = RNA_PATH + 'negative_numbers.pdb'

# RNA with a single large strand break
OPPOSITEGAP =TEST_DATA_PATH+'rna_structures/rna_huge_gap.pdb'
RNA_TWO_PIECES = TEST_DATA_PATH+'rna_structures/rna_2pieces.ent'

# RNA with a different kind of backbone break in each residue.
BROKEN_BACKBONE = TEST_DATA_PATH+'rna_structures/broken_backbones.ent'
BB_MESSED_UP = TEST_DATA_PATH+'rna_structures/bb_messed_up.pdb'

# RNA with a clash
CLASHING_STRUCTURE = TEST_DATA_PATH+'clash/clash.ent'

# RNA with missing atoms
INCOMPLETE_BACKBONE = TEST_DATA_PATH+'rna_structures/incomplete_backbone.pdb'

# fragments
FRAGMENT1 = TEST_DATA_PATH+'rna_structures/fragment1.pdb'
FRAGMENT2 = TEST_DATA_PATH+'rna_structures/fragment2.pdb'
SMALL_FRAGMENT = TEST_DATA_PATH+'rna_structures/small_fragment.pdb'
BROKEN_FRAGMENT = TEST_DATA_PATH+'rna_structures/broken_fragment.pdb'
UNKNOWN_FRAGMENT = TEST_DATA_PATH+'rna_structures/unknown_fragment.pdb'
FCCD_EXAMPLE = TEST_DATA_PATH+'rna_structures/fccd_example.ent'
FIXABLE_BACKBONE =TEST_DATA_PATH+'rna_structures/bb_to_fix.pdb'
BROKEN_HELIX_EXAMPLE = TEST_DATA_PATH+'rna_structures/broken_helix.pdb'
BAD_LOOP_INSERTION = TEST_DATA_PATH+'rna_structures/bad_loop_insertion.pdb'
MISSING_PHOSPHATES = TEST_DATA_PATH+'rna_structures/missing_phosphates.pdb'
MISSING_PHOSPHATES2 = TEST_DATA_PATH+'rna_structures/missing_p_17.pdb'

#
# other structure examples
#

# SO4 is not a modification
SULFATE = TEST_DATA_PATH+'other/sulfate.pdb'

# DNA strand, chain E
DNA_WITH_MISMATCH = TEST_DATA_PATH+'other/dna_with_mismatch.pdb'

# Protein is not RNA
PROTEIN_STRUCTURE = TEST_DATA_PATH+'other/protein.pdb'

# nasty PDB structure for check_pdb tests
NASTY_PDB = TEST_DATA_PATH+'check_pdb/1EHZ2.pdb'
PDB_WITH_LIGAND = TEST_DATA_PATH+'check_pdb/3FU2.pdb'
PDB_WITH_DOUBLE_COORDINATES = TEST_DATA_PATH+'check_pdb/1L3D.pdb'
PDB_WITH_AMP = TEST_DATA_PATH+'check_pdb/1RAW.pdb'
PDB_WITH_PHOSPHATE = PDB_WITH_LIGAND
PDB_WITHOUT_DOUBLE_COORD = TEST_DATA_PATH+'check_pdb/disordered_residues.pdb'
# bulge motif
BULGE_MOTIF = TEST_DATA_PATH+ 'rna_structures/bulge_motif.pdb'

#
# RNA with unknown residues to test mutations (manualy modifierd fragment of 1ehz, 12 resi)
#
PDB_UNK = TEST_DATA_PATH+ 'rna_structures/unk.pdb'


#
# strange PDBs for testing PdbConverter (*,HETATM itd.)
#
# example from JMB
PDB_STRANGE = TEST_DATA_PATH+'rna_structures/hybrid2.pdb'

BAD_TEMPLATE = TEST_DATA_PATH+'geometry/mini_bad_geometry.pdb'
