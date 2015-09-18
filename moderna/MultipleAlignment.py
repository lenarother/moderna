#!/usr/bin/env python

"""
Multiple-alignment based modeling for ModeRNA
"""

from sequence.RNAAlignment import RNAAlignment
from sequence.ModernaSequence import Sequence
from Template import Template
from RNAModel import RnaModel
from moderna import create_model
from math import sqrt, sin, acos
import re
import Bio


DISTANCE_THRESHOLD = 3.0    # C4' atom distance threshold (angstroms)

def get_atom_distance(coords1, coords2):
    """
    Returns the distance between two atoms in 3D space.
    """
    return sqrt(pow(coords1[0] - coords2[0], 2) +\
                pow(coords1[1] - coords2[1], 2) +\
                pow(coords1[2] - coords2[2], 2))

def centroid(*coords):
    """
    Returns centroid (x,y,z) for the specified coordinate set.
    It takes a list of 3D coordinates as an input.
    """
    x = y = z = 0
    for c in coords:
        x+=c[0]
        y+=c[1]
        z+=c[2]
    N = len(coords)
    x = x / N
    y = y / N
    z = z / N
    return x, y, z

def parse_alignment(seq_alignment):
    """
    Returns all sequences from the alignment as a tuple.
    """
    spl_seq = seq_alignment.split(">")
    spl_seq.pop(0)
    return [(("").join(seq.split("\n")[0]), ("").join(seq.split("\n")[1:]))\
            for seq in spl_seq]


class IncorrectAlignmentError(Exception):
    pass


class ErroneousStructure(Exception):
    pass


class NoTemplateError(Exception):
    pass


class MultipleAlignment(object):
    """
    Stores the multiple alignment of RNA secondary structures
    and sequences together with accompanying templates.
    """
    
    def __init__(self, seq_alignment, templates):
        self.seq_alignment = seq_alignment
        self.seq_alignment_tuple = parse_alignment(seq_alignment)
        self.template_matrix = []
        self.set_templates(templates)
    
    def __repr__(self):
        return "<RNA secondary structure alignment>\n" +\
               "\n".join(self.sec_struct_alignment)
    
    def all_template_numbers(self):
        """
        Returns all template numbers in the alignment.
        """
        return range(self.temp_count)
    
    def get_pairs(self):
        """
        Returns pairs for each residue from all templates.
        -1 means no pairing. None means gap.
        """
        return self.pairs

    def get_secstruct(self, template, seq_with_gaps=None):
        """
        Returns secondary structure of a template.
        """
        import os, time
        combined = []

        suffix = time.time()
        X3DNA_OUT = "x3dna.%s.out" % suffix
        tmp_file = "x3dna.%s.pdb" % suffix

        gap_numbers = []
        if seq_with_gaps is not None:
            numbers = []
            length = len(seq_with_gaps)
            for i in range(length):
                if seq_with_gaps[i]!="-":
                    numbers.append(i)
                else:
                    gap_numbers.append(i)
        else:
            length = len(template)
            numbers = range(length)

        m = RnaModel()
        i = 0
        for r in template:
            m.add_residue(r, number=str(numbers[i]))
            i += 1
        m.write_pdb_file(tmp_file)

        reg_exp = re.compile('>A:[.]*([0-9]+)_:[^:]*:[.]*([0-9]+)_:A<')

        exec_path = os.path.join(os.path.dirname(__file__), "x3dna.sh")
        os.system("%s %s %s" % (exec_path, tmp_file, X3DNA_OUT))

        try:
            fh = open(X3DNA_OUT).read()
            pairs = reg_exp.findall(fh)
            pairs = [(int(p[0]), int(p[1])) for p in pairs]
        except:
            pairs = []

        fnames = ['bp_order.dat', 'col_chains.scr', 'col_helices.scr', 'ref_frames.dat',
                  'bestpairs.pdb', 'hel_regions.pdb', tmp_file, X3DNA_OUT]

        for fn in fnames:
            if os.path.exists(fn):
                os.remove(fn)

        vienna_tab = ["."] * length
        for i in gap_numbers:
            vienna_tab[i] = "-"
        for pair in pairs:
            vienna_tab[min(pair)] = "("
            vienna_tab[max(pair)] = ")"
        vienna = "".join(vienna_tab)

        return {"pairs": pairs, "gaps": gap_numbers, "vienna": vienna}
    
    def get_pairing(self, vienna, pairs=None, gaps=None):
        """
        Returns a list with the positions of paired nucleotides,
        -1 if a nucleotide has no pair, or None for a gap.
        """
        stack = []
        length = len(vienna)
        pair_table = [-1] * length
        if pairs is not None and gaps is not None:
            for gap in gaps:
                pair_table[gap] = None
            for pair in pairs:
                res1 = min(pair)
                res2 = max(pair)
                pair_table[res1] = res2
                pair_table[res2] = res1
        else:
            for i in range(length):
                if vienna[i]=="(":
                    stack.append(i)
                elif vienna[i]==")":
                    try:
                        prev_pos = stack.pop()
                    except IndexError:
                        raise ErroneousStructure
                    pair_table[prev_pos] = i
                    pair_table[i] = prev_pos
                elif vienna[i]=="-":
                    pair_table[i] = None
        return pair_table

    def make_sec_struct_alignment(self):
        """
        Creates the secondary structure alignment of the templates.
        Updates internal variables with the new data.
        """

        templates = self.get_templates()
        seq_alignment = self.get_seq_alignment_tuple()
        if len(seq_alignment)>0:
            seq_aln_len = len(seq_alignment[0][1])
        else:
            seq_aln_len = 0
        
        ss_alignment = []
        sec_structs = []
        for i, template in enumerate(templates):
            sec_struct = self.get_secstruct(template, seq_alignment[i+1][1])
            sec_structs.append(sec_struct)
            ss_alignment.append(sec_struct["vienna"])
        self.sec_struct_alignment = ss_alignment
        
        self.temp_count = len(self.sec_struct_alignment)
        if self.temp_count>0:
            self.length = len(self.sec_struct_alignment[0])
        else:
            self.length = 0
        self.pairs = []
        for ss in sec_structs:
            self.pairs.append(self.get_pairing(**ss))

        self.guides = self.get_guides()

    def get_sec_struct_alignment(self):
        """
        Returns secondary structure alignment.
        """
        return self.sec_struct_alignment
   
    def get_seq_alignment(self):
        """
        Returns sequence alignment.
        """
        return self.seq_alignment

    def get_seq_alignment_tuple(self):
        """
        Returns sequence alignment as a tuple.
        """
        return self.seq_alignment_tuple
    
    def set_templates(self, template_list):
        """
        Sets the template list for the alignment.
        Templates are associated with respective structures,
        i.e. the 3rd structure in the alignment describes the 3rd template.
        """
        self.templates = template_list
        self.make_sec_struct_alignment()

    def get_templates(self):
        """
        Returns all templates.
        """
        return self.templates
    
    def is_correct_alignment(self, alignment):
        """
        Checks the length (all structures should have the same)
        and presence of illegal characters in the alignment.
        """
        length = None
        for line in alignment:
            if length is None:
                length = len(line)
            else:
                if len(line)!=length:
                    return False
            if re.match(r"[^\(\)\.\-]", line):
                return False
        return True
    
    def get_guides(self):
        """
        Returns a list with "guide" templates for each position.
        """
        guides = [None] * self.length
        for pos in range(self.length):
            for temp_num in self.all_template_numbers():
                if self.sec_struct_alignment[temp_num][pos]!="-":
                    guides[pos] = temp_num
                    break
        return guides
        
    def get_structure_matrix(self):
        """
        Returns a binary matrix, denoting which templates shall be used
        at specified position, basing on their secondary structures.
        """
        matrix = []
        for i in self.all_template_numbers():
            matrix.extend([[0] * self.length])
        for temp_num in self.all_template_numbers():
            for pos in range(self.length):
                try:
                    if self.pairs[temp_num][pos] is not None and\
                       self.pairs[temp_num][pos] in\
                       (-1, self.pairs[self.guides[pos]][pos]):
                        matrix[temp_num][pos] = 1
                except TypeError: # no guide at that position
                    matrix[temp_num][pos] = 0
        return matrix

    def get_atoms(self, temp_num, name):
        """
        Returns the coordinates of all atoms with the specified name
        in the template. If two names are given separated by slash ("/"),
        the second one is used when the residue is a pyrimidine.
        """
        name2 = None
        if "/" in name:
            names = name.split("/")
            name = names[0]
            name2 = names[1]
        template = self.templates[temp_num]
        atoms = []
        for residue in template:
            if residue.pyrimidine and name2 in residue:
                atoms.append(residue[name2])
            elif name in residue:
                atoms.append(residue[name])
        coords = [atom.get_coord() for atom in atoms]
        for pos in range(len(self.sec_struct_alignment[temp_num])):
            if self.sec_struct_alignment[temp_num][pos]=="-":
                coords.insert(pos, None)
        return coords
    
    def get_template_matrix(self):
        """
        Returns binary matrix, denoting which templates shall be used
        at specified position, basing on their secondary structures
        AND the distances of C4' atoms from the guide template.
        """     
        if len(self.template_matrix)>0:
            return self.template_matrix
        struct_matrix = self.get_structure_matrix()
        dist_matrix = []
        self.coords = {"C4'": [], "C1'": [], "C2'": [],
                       "N9/N1": [], "C8/C6": [], "C4/C2": []}
        try:
            for name in self.coords.keys():
                self.coords[name] = [self.get_atoms(temp_num, name) for\
                                     temp_num in self.all_template_numbers()]
        except IndexError:
            raise NoTemplateError
        guide_coords = []
        for pos in range(self.length):
            try:
                guide_coords.append(self.coords["C4'"][self.guides[pos]][pos])
            except TypeError: # no guide at that position
                guide_coords.append(None)
        self.dists = [] #XXX
        for temp_num in self.all_template_numbers():
            tmp = [0] * self.length
            for pos in range(self.length):
                distance = None
                if self.guides[pos]==temp_num:
                    tmp[pos] = 1
                elif struct_matrix[temp_num][pos]==1:
                    distance = get_atom_distance(self.coords["C4'"][temp_num][pos],\
                               guide_coords[pos])
                    if distance<=DISTANCE_THRESHOLD:
                        tmp[pos] = 1
                if distance is None: distance = "-" #XXX
                if temp_num==1: self.dists.append(str(distance)) #XXX
            dist_matrix.append(tmp)
        self.template_matrix = dist_matrix
        return dist_matrix

    def target_atom_coords(self, name):
        """
        Returns a list of coordinates for atoms with the specified name
        for each nucleotide that can be resolved from the templates.
        """
        matrix = self.get_template_matrix()
        coords = []
        for pos in range(self.length):
            selected = []
            for temp_num in self.all_template_numbers():
                if matrix[temp_num][pos]==1:
                    selected.append(self.coords[name][temp_num][pos])
            if len(selected)>0:
                coord = centroid(*selected)
            else:
                coord = None
            coords.append(coord)
        return coords


class AlignmentWrapper(RNAAlignment):
    """
    Wrapper for RNAAlignment class, allowing for using multiple alignments
    in ModeRNA.
    Input: MultipleAlignment object.
    """
    
    def __init__(self, multaln, template_seq=None):
        seq_alignment_tab = multaln.get_seq_alignment_tuple()
        single_alignment = self.shrink_to_single_alignment(seq_alignment_tab,\
                           multaln.get_guides(), template_seq)
        ssa = parse_alignment(single_alignment)
        RNAAlignment.__init__(self, ssa[0][0], Sequence(ssa[0][1]),
                                    ssa[1][0], Sequence(ssa[1][1]))
   
    def shrink_to_single_alignment(self, alignment, guides, template_seq=None):
        """
        "Shrinks" the alignment to a single target-template alignment
        with template sequence either specified explicitly (template_seq)
        or created from multiple alignment's guide list.
        """
        if template_seq:
            temp_seq = template_seq
        else:
            temp_seq = ""
            length = len(alignment[0][1])
            for pos in range(length):
                if guides[pos] is not None:
                    temp_seq += alignment[guides[pos]+1][1][pos]
                else:
                    temp_seq += "-"
        single_alignment = (">" + alignment[0][0], alignment[0][1],\
                            ">template", temp_seq)
        return ("\n").join(single_alignment)


class TemplateWrapper(Template):
    """
    Wrapper for Template class, allowing for using multiple templates
    in ModeRNA.
    Input: MultipleAlignment object
    """
    
    def __init__(self, multaln):
        templates = multaln.get_templates()
        struct_matrix = multaln.get_structure_matrix()
        sec_structs = multaln.get_sec_struct_alignment()
        temp_list = []
        for temp_num in range(len(templates)):
            iterate = templates[temp_num].__iter__()
            tmp = []
            for pos in range(len(struct_matrix[temp_num])):
                if sec_structs[temp_num][pos]!="-":
                    res = iterate.next()
                if struct_matrix[temp_num][pos]==1:
                    tmp.append(res)
                else:
                    tmp.append(None)
            temp_list.append(tmp)
        template = RnaModel()
        target_seq = multaln.get_seq_alignment_tuple()[0][1]
        target_purines = [base in ("A", "G") for base in target_seq]
        self.guides = multaln.get_guides()
        self.real_sequence = ""
        list_of_used_templates = []
        pair_list = []
        for pos in range(len(self.guides)):
            if self.guides[pos] is not None:
                temp_used = None
                if temp_list[self.guides[pos]][pos].purine==target_purines[pos]:
                    temp_used = self.guides[pos]
                if temp_used is None:
                    for temp_num in range(len(temp_list)):
                        if temp_list[temp_num][pos] is not None\
                        and temp_list[temp_num][pos].purine==target_purines[pos]:
                            temp_used = temp_num
                            break
                if temp_used is None:
                    temp_used = self.guides[pos]
                # if this is a "closing bracket",
                # use the same template as with matching one:
                opening_pos = self.get_opening_bracket(pos, multaln, temp_used)
                if opening_pos is not None:
                    which_template_there = list_of_used_templates[opening_pos]
                    if temp_list[which_template_there][pos] is not None:
                        # avoid setting temp_used to None
                        temp_used = which_template_there
                        pair_list.append((opening_pos, pos))
                list_of_used_templates.append(temp_used)
                res = temp_list[temp_used][pos].copy()
                res.change_number(str(pos))
                template.add_residue(res)
                self.real_sequence += res.original_base
            else:
                list_of_used_templates.append(None)
                self.real_sequence += "-"
        Template.__init__(self, template.get_structure(),
                          template_data_type="structure",
                          template_chain_name=template.chain_name)
        back_coords = {"C4'": multaln.target_atom_coords("C4'"),
                       "C1'": multaln.target_atom_coords("C1'"),
                       "C2'": multaln.target_atom_coords("C2'"),}
        self.replace_sugars(back_coords, temp_list)
        self.move_residues(back_coords, pair_list)
        base_coords = {"N9/N1": multaln.target_atom_coords("N9/N1"),
                       "C8/C6": multaln.target_atom_coords("C8/C6"),
                       "C4/C2": multaln.target_atom_coords("C4/C2"),}
        self.rotate_bases(base_coords)
        for pos in range(len(struct_matrix[0])): #XXX
            r1 = temp_list[0][pos]
            if r1 is None:
                id1 = "-"
                name1 = "-"
            else:
                id1 = str(r1.id[1])
                name1 = r1.resname
            r2 = temp_list[1][pos]
            if r2 is None:
                id2 = "-"
                name2 = "-"
            else:
                id2 = str(r2.id[1])
                name2 = r2.resname
            print "\t".join((id1, name1, id2, name2, multaln.dists[pos]))

    def get_opening_bracket(self, position, multaln, temp_used):
        """
        Gets opening bracket (matching pair) for given position
        or returns None if there's no closing bracket at that position.
        """
        if temp_used is None:
            return None
        ss = multaln.get_sec_struct_alignment()
        if ss[temp_used][position]!=")":
            return None
        pairs = multaln.get_pairs()
        return pairs[temp_used][position]
    
    def replace_sugars(self, centroids, repl_templates):
        """
        Replaces sugars in all residues with the sugars
        closest to given centroid coordinates.
        """
        num_templates = len(repl_templates)
        i = 0
        for r in self:
            while self.real_sequence[i]=="-":
                i += 1
            d = []
            for t in repl_templates:
                if t[i] is not None:
                    dist1 = get_atom_distance(t[i]["C4'"].coord,\
                                              centroids["C4'"][i])
                    dist2 = get_atom_distance(t[i]["C1'"].coord,\
                                              centroids["C1'"][i])
                    dist3 = get_atom_distance(t[i]["C2'"].coord,\
                                              centroids["C2'"][i])
                    d.append(sqrt(dist1 + dist2 + dist3))
                else:
                    d.append(1e4)
            repl_res = repl_templates[d.index(min(d))][i]
            for atom in repl_res:
                if not "'" in atom.name: continue
                r[atom.name].set_coord(atom.coord)
            i += 1
    
    def rotate_bases(self, new_coords):
        """
        Rotates each base superimposing it against (N9, C8, C4) centroids
        for purines and (N1, C6, C2) for pyrimidines.
        """
        i = 0
        exclude_atoms = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'",\
                         "C3'", "O3'", "C2'", "O2'", "C1'")
        for r in self:
            while self.guides[i] is None:
                i += 1
            centroids = []
            for name in ("N9/N1", "C8/C6", "C4/C2"):
                centroids.append(Bio.PDB.Atom.Atom(name, new_coords[name][i],\
                                 0.0, 1.0, " ", " %s " % name, i, name[0]))
            sup = Bio.PDB.Superimposer()
            if r.pyrimidine:
                sup.set_atoms(centroids, (r["N1"], r["C6"], r["C2"]))
            else:
                sup.set_atoms(centroids, (r["N9"], r["C8"], r["C4"]))
            rot, tran = sup.rotran
            for atom in r:
                if not atom.name in exclude_atoms:
                    atom.transform(rot, tran)
            i += 1
    
    def move_residues(self, new_coords, pair_list):
        """
        Superimposes all residues against (C4', C1', C2') centroids.
        Paired residues are treated as a whole.
        """
        i = 0
        openings, closings = zip(*pair_list)
        transform_data = {} # contains "opening bracket" residues and their rot/tran data
        for r in self:
            while self.guides[i] is None:
                i += 1
            centroids = []
            for name in ("C4'", "C1'", "C2'"):
                centroids.append(Bio.PDB.Atom.Atom(name, new_coords[name][i],\
                                 0.0, 1.0, " ", " %s " % name, i, name[0]))
            sup = Bio.PDB.Superimposer()
            sup.set_atoms(centroids, (r["C4'"], r["C1'"], r["C2'"]))
            rot, tran = sup.rotran
            if i in openings:
                transform_data[i] = (r, rot, tran)
            elif i in closings:
                index = closings.index(i)
                r_op, rot_op, tran_op = transform_data[openings[index]]
                rot = (rot + rot_op) / 2
                tran = (tran + tran_op) / 2
                r_op.transform(rot, tran)
            if not i in openings:
                r.transform(rot, tran)
            i += 1

