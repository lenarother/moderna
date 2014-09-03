#!/usr/bin/env python
#
# test_geometry_statistics.py
#
# unit tests for checkinh geometry statisticcs
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Kristian Rother"
__copyright__ = "Copyright 2008, Genesilico"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "alpha"

from unittest import main, TestCase
from moderna.analyze.GeometryStatistics import GeometryStatistics, PDBSetGeometryStatistics, GeometryResult, GeometryExpression, AtomDefinition
from moderna.ModernaStructure import ModernaStructure
from moderna.RNAResidue import RNAResidue
from test_data import *
import os

EXPRESSIONS  = (
    "X:C5',X:C4'",
    "X:C5',X:C4',X:C3'",
    "X:C5',X:C4',X:C3',X:O3'",
    "NUC:X[G]:N1,X[A]:N1,Y:O2'",
    "NUC:Y[C]:O2',NUC:Y+1[C]:O2'",
)

model = ModernaStructure('file', RNA_1EHZ)

residue_list = [r for r in model]
"""
<Residue   G het=  resseq=1 icode= >
<Residue   C het=  resseq=2 icode= >
<Residue   G het=  resseq=3 icode= >
<Residue   G het=  resseq=4 icode= >
"""

class AtomDefinitionTests(TestCase):
    
    def setUp(self):
        self.plain = AtomDefinition("X:C5'")
        self.resname = AtomDefinition("X[C]:N1")
        self.chain = AtomDefinition("X{A}:N1")
        self.chain_combo = AtomDefinition("X+1{A}[G]:N1")
        self.bad_chain = AtomDefinition("X{B}:N1")
        self.chaintype = AtomDefinition("NUC:X:O2'")
        self.offset = AtomDefinition("X+1:C5'")
        self.neg_offset = AtomDefinition("X-1:C5'")
        self.combined = AtomDefinition("NUC:Y+1[G]:O2'")
        
        self.purine = AtomDefinition("X[R]:N1")
        self.pyrimidine = AtomDefinition("X[Y]:N1")
        
    def test_get_atom(self):
        atom = self.plain.get_atom(residue_list,2)
        self.assertEqual(atom.name,"C5'")
        self.assertEqual(atom.parent.id[1],3)
        self.assertEqual(atom.parent.resname,'G')
        # check whether the +1 residue is identified
        atom = self.offset.get_atom(residue_list,2)
        self.assertEqual(atom.parent.id[1],4)
        atom = self.offset.get_atom(residue_list,len(residue_list)-1)
        self.assertEqual(atom,None)
        # check whether the -1 residue is identified
        atom = self.neg_offset.get_atom(residue_list,2)
        self.assertEqual(atom.parent.id[1],2)
        atom = self.neg_offset.get_atom(residue_list,0)
        self.assertEqual(atom,None)
        # check whether the residue name is taken into account
        atom = self.resname.get_atom(residue_list,1)
        self.assertEqual(atom.parent.resname,'C')
        atom = self.resname.get_atom(residue_list,2)
        self.assertEqual(atom,None)
        # check whether the chain identifiers are used
        # NO CHAIN IDENTIFIERS IN MODERNA
        #atom = self.chain.get_atom(residue_list,1)
        #self.assertEqual(atom.parent.resname,'C')
        #atom = self.chain_combo.get_atom(residue_list,1)
        #self.assertEqual(atom.parent.resname,'G')
        #atom = self.bad_chain.get_atom(residue_list,1)
        #self.assertEqual(atom,None)

    def test_purine(self):
        atom = self.purine.get_atom(residue_list,0)
        self.assertNotEqual(atom,None)
        atom = self.purine.get_atom(residue_list,1)
        self.assertEqual(atom,None)

    def test_pyrimidine(self):
        atom = self.pyrimidine.get_atom(residue_list,1)
        self.assertNotEqual(atom,None)
        atom = self.pyrimidine.get_atom(residue_list,2)
        self.assertEqual(atom,None)
        
        
        
class GeometryExpressionTests(TestCase):

    def setUp(self):
        self.ge = GeometryExpression("X:P,Y:C4'")
        
    def test_init(self):
        """Checks whether the right number of AtomDefinition objects is found."""
        for ex,size in zip(EXPRESSIONS,[2,3,4,3,2]):
            ge = GeometryExpression(ex)
            self.assertEqual(len(ge),size)
            
    def test_get_residue_codes(self):
        """The residue variables should be extracted correctly."""
        for ex, codes in zip(EXPRESSIONS,[['X'],['X'],['X'],['X','Y'],['Y']]):
            ge = GeometryExpression(ex)
            self.assertEqual(ge.get_residue_codes(),codes)
            
            
    def test_get_residues_by_codes(self):
        """Should return lists of residues for each code."""
        residues = self.ge.get_residues_by_codes(model)
        self.assertEqual(len(residues),2)
        for r in residues[0]+residues[1]:
            self.assertTrue(isinstance(r, RNAResidue))
        
    def test_has_distinct_residues(self):
        """Checks for non-duplicate indices in the list."""
        residues = (('a','b','c','d'),('1','2','a','b','3'))
        self.assertTrue(self.ge.has_distinct_residues(residues,[0,0]))
        self.assertTrue(self.ge.has_distinct_residues(residues,[1,4]))
        self.assertTrue(self.ge.has_distinct_residues(residues,[3,3]))
        self.assertTrue(self.ge.has_distinct_residues(residues,[1,2]))
        self.assertFalse(self.ge.has_distinct_residues(residues,[0,2]))
        self.assertFalse(self.ge.has_distinct_residues(residues,[1,3]))
        
    def test_get_atom_combination(self):
        """Should return the indexed atoms"""
        residues = (residue_list[:20],residue_list[10:30]) 
        examples = (
            ((0,0),(residue_list[0],residue_list[10])),
            ((2,1),(residue_list[2],residue_list[11])),
            ((3,2),(residue_list[3],residue_list[12])),
            )
        for r_index, resis in examples:
            atoms = self.ge.get_atom_combination(residues,r_index)
            matches = (resis[0]["P"],resis[1]["C4'"])
            for i in range(2):
                self.assertEqual(atoms[i],matches[i])
        
    def test_check_atom_combination(self):
        """Should check whether the atoms are matched by AtomDefinitions."""
        self.assertTrue(self.ge.check_atom_combination([1,2,3,4]))
        self.assertFalse(self.ge.check_atom_combination([1,2,None,4]))
            
    def test_get_atoms(self):
        i = 0
        for atoms in self.ge.get_atoms(model):
            self.assertEqual(len(atoms),2)
            i += 1
        self.assertTrue(i>0)
        
        
    
class GeometryStatisticsTests(TestCase):
    def setUp(self):
        t = ModernaStructure('file', MINI_TEMPLATE)
        self.gs = GeometryStatistics(t)
        
    def test_main_part(self):
        """GeometryStatistics should calculate dist, angle and dihedral data."""
        result = self.gs.get_distances("X:C5',X:C4'")
        self.assertTrue(len(result)>0)
        result = self.gs.get_angles("X:C5',X:C4',X:C3'")
        self.assertTrue(len(result)>0)
        result = self.gs.get_dihedrals("X:C5',X:C4',X:C3',X:O3'")
        self.assertTrue(len(result)>0)
        
    def test_get_distances(self):
        result = self.gs.get_distances("X:C5',X:C4'")
        self.assertTrue(result)
    
    def test_get_angles(self):
        result = self.gs.get_angles("X:C5',X:C4',X:C3'")
        self.assertTrue(result)
    
    def test_get_dihedrals(self):
        result = self.gs.get_dihedrals("X:C5',X:C4',X:C3',X:O3'")
        self.assertTrue(result)

class PDBSetGeometryStatisticsTests(TestCase):
    def setUp(self):
        self.gs = PDBSetGeometryStatistics(TEST_DATA_PATH+'geometry')
        if os.access('table.txt',os.F_OK): os.remove('table.txt')
        if os.access('plot.png',os.F_OK): os.remove('plot.png')

    def tearDown(self):
        if os.access('table.txt',os.F_OK): os.remove('table.txt')
        if os.access('plot.png',os.F_OK): os.remove('plot.png')

    def test_main_part(self):
        """GeometryStatistics should calculate dist, angle and dihedral data."""
        result = self.gs.get_distances("X:C5',X:C4'")
        self.assertTrue(len(result)>0)
        result = self.gs.get_angles("X:C5',X:C4',X:C3'")
        self.assertTrue(len(result)>0)
        result = self.gs.get_dihedrals("X:C5',X:C4',X:C3',X:O3'")
        self.assertTrue(len(result)>0)
        
    def test_get_structures(self):
        i = 0
        for s in self.gs.get_structures():
            self.assertTrue(isinstance(s,ModernaStructure))
            i += 1
        self.assertEqual(i,2)
        
    def test_tabular_output(self):
        """Results should be written to a text file."""
        result = self.gs.get_distances("X:C5',X:C4'")
        result.write_table('table.txt')
        self.assertTrue(os.access('table.txt',os.F_OK))

    def _test_plotting(self):
        """NOT IMPLEMENTED:Results should be written to a png plot."""
        result = self.gs.get_distances("X:C5',X:C4'")
        result.write_plot('plot.png')
        self.assertTrue(os.access('plot.png',os.F_OK))
        
    def test_get_distances(self):
        result = self.gs.get_distances("X:C5',X:C4'")
        self.assertTrue(result)
    
    def test_get_angles(self):
        result = self.gs.get_angles("X:C5',X:C4',X:C3'")
        self.assertTrue(result)
    
    def test_get_dihedrals(self):
        result = self.gs.get_dihedrals("X:C5',X:C4',X:C3',X:O3'")
        self.assertTrue(result)
    
        

class GeometryResultTests(TestCase):
    def setUp(self):
        self.plain = GeometryResult('plain data')
        self.angles = GeometryResult('angles', angles=True)
        [self.plain.append(x) for x in [(2.0,'A1', 'B1'), (2.4,'A1', 'B1'), (2.1,'A1', 'B1'), (1.7,'A1', 'B1'), (1.8,'A1', 'B1')]]
        [self.angles.append(x) for x in [(90.0,'A1', 'B1', 'C2'), (350.0,'A1', 'B1', 'C2'), (270.0,'A1', 'B1', 'C2'), (50.0,'A1', 'B1', 'C2')]]
        if os.access('table.txt',os.F_OK): os.remove('table.txt')
        if os.access('plot.png',os.F_OK): os.remove('plot.png')

    def tearDown(self):
        if os.access('table.txt',os.F_OK): os.remove('table.txt')
        if os.access('plot.png',os.F_OK): os.remove('plot.png')

    def test_length(self):
        self.assertEqual(len(self.plain),5)
        self.assertEqual(len(self.angles),4)
        
    def test_avg(self):
        self.assertAlmostEqual(self.plain.get_avg(), 2.0)
        self.assertAlmostEqual(self.angles.get_avg(),10.0)
        
    def test_stddev(self):
        self.assertAlmostEqual(self.plain.get_stddev(), 0.273861278)
        self.assertAlmostEqual(self.angles.get_stddev(), 78.315600829)
        
    def test_tabular_output(self):
        self.plain.write_table('table.txt')
        self.assertTrue(os.access('table.txt',os.F_OK))    

    def _test_plotting(self):
        """NOT IMPLEMENTED"""
        self.angles.write_plot('plot.png')
        self.assertTrue(os.access('plot.png',os.F_OK))
    


if __name__ == '__main__':
    main()
    
