================
ModeRNA Commands
================

These are usage examples of ModeRNA functions.


add_modification
================

.. code-block:: python

    # add a 2-methyl-Guanosine in position 63.
    m  = load_model('1QF6_B_tRNA.pdb', 'B')
    add_modification(m['63'], 'm2G')

add_all_modifications
=====================

.. code-block:: python

    add_all_modifications(t,a,m)


add_pair_to_base
================

.. code-block:: python

    add_pair_to_base(m, '2', '20', 'C')


analyze_geometry
================

.. code-block:: python
    
    analyze_geometry(t)
    analyze_geometry(m)

analyze_geometry
================

.. code-block:: python

    t = load_template('1QF6_B_tRNA.pdb', 'B')
    analyze_geometry(t)

apply_indel
===========

.. code-block:: python

    # prepare a model
    t = load_template('1QF6_B_tRNA.pdb', 'B')
    m = create_model()
    copy_some_residues(t['31':'35']+t['38':'42'],m)

    # Find the best fitting fragment for the missing two
    apply_indel(m, '35', '38', 'CA')''',


apply_alignment
===============

.. code-block:: python

    apply_alignment(t,a,m)

apply_all_indels
================

.. code-block:: python

    apply_all_indels(a,m)

apply_missing_ends
==================

.. code-block:: python

    apply_missing_ends(a, m)

change_sequence
===============

.. code-block:: python

    change_sequence(m, 'AGCUAGCU')

clean_structure
===============

.. code-block:: python

    # cleaning up a loaded template:
    # removes water, ions, amino acids, and unknown residues
    # replaces O1P and O2P in atom names by OP1 and OP2
    # replaces * in atom names by '
    clean_structure(t)

create_model
============

.. code-block:: python

    m = create_model()
    m = create_model(model_chain_name='K')
    m = create_model(t,a)
    m = create_model(t,a,model_chain_name='K')

copy_identical_residues
=======================

.. code-block:: python

    copy_identical_residues(t,a,m)
    copy_identical_residues(t,a,m,strict=0)
    copy_identical_residues(t,a,m,modifications=0)

copy_single_residue
===================

.. code-block:: python

    copy_single_residue(t['3'],m)
    copy_single_residue(t['5A'],m)
    copy_single_residue(t['3'],m,'15')
    copy_single_residue(m['5A'],m,'100')
    copy_single_residue(t['3'],m,'3B',strict=0)

copy_some_residues
==================

.. code-block:: python

    copy_some_residue(t['3':'5'],m)
    copy_some_residue(t['3':'5A'],m)
    copy_some_residue([t['3'],t['7'],t['8']],m)
    copy_some_residue([t['3'],t['7']],m,new_numbers=['100','101'])
    copy_some_residue(t['3':'12'],m,strict=0)

create_fragment
===============

.. code-block:: python

    f = create_fragment('single_strand.pdb',anchor5=m['20'])
    f = create_fragment('single_strand.pdb', chain_name='A', anchor3=m['1'])
    f = create_fragment('hairpin.pdb', anchor5=m['12'], anchor3=m['15'])
    f = create_fragment('hairpin.pdb', anchor5=m['12'], anchor3=m['15'], sequence='AG')

delete_residue
==============

.. code-block:: python

    delete_residue('6',m)

examine_structure
=================

.. code-block:: python

    # examine a loaded template for irregularities:
    examine_structure(t)
    examine_structure(t,'logfile.log')


exchange_single_base
====================

.. code-block:: python

    exchange_single_base(m['3'],'C')
    exchange_single_base(t['3'],'G',m)
    exchange_single_base(t['3'],'G',m,new_number='5A')
    exchange_single_base(t['3'],'G') # modifies the template!


exchange_some_bases
===================

.. code-block:: python

    exchange_some_bases(m['3':'5'],['C','U','G'])
    exchange_some_bases(t['3':'5'],['C','U','G'],m)
    exchange_some_bases([t['3'],t['10']],['A','G'],m,['103','106'])
    exchange_some_bases(t['3':'5'],['C','U','G']) # modifies the template!

exchange_mismatches
===================

.. code-block:: python

    exchange_mismatches(t,a,m)

extend_helix
============

.. code-block:: python

    extend_helix(m, '2', '37', 'CGAA_UUCG')

find_fragment
=============

.. code-block:: python

    candidates = find_fragment(m,'7','12','AGGU', 10)
    insert_fragment(m, candidates[0])

find_modifications
==================

.. code-block:: python

    find_modifications(t)
    mods = find_modifications(m)


find_clashes
============

.. code-block:: python

    find_clashes(m)
    find_clashes([m['5'], m['6'], m['7']])
    pairs = find_clashes(m)

find_fragment
=============

.. code-block:: python

    candidates = find_fragment(m, '7', '12', 'AGCU', 20)
    insert_fragment(m, candidates[0])

find_modifications
==================

.. code-block:: python

    # find modifications in any structure file
    t = load_template('1QF6_B_tRNA.pdb', 'B')
    print find_modifications(t)


fix_backbone
============

.. code-block:: python

    m = load_model('broken_bb.pdb','A')
    print m.get_sequence()

    # check and fix the entire model
    fix_backbone(m)
    print m.get_sequence()

    # check and fix the connection between residues 4 and 5
    fix_backbone(m, '4', '5')


get_base_pairs
==============

.. code-block:: python

    bp = get_base_pairs(struc)
    print bp

get_sequence
============

.. code-block:: python

    get_sequence(t)
    seq = get_sequence(m)


get_secstruc
============

.. code-block:: python

    get_secstruc(t)
    ss = get_secstruc(m)

get_stacking
============

.. code-block:: python

    stacking = get_stacking(struc)
    print stacking

insert_fragment
===============

.. code-block:: python

    insert_fragment(m,f)

insert_two_strand_fragment
==========================

.. code-block:: python

    insert_two_strand_fragment(m, '2', '37', '5', '34',\
    '101', '120', '107', '114', 'my_fragment.pdb', 'A')
    # fragment candidates returned by find_fragment go as well
    insert_fragment(m, candidates[2])

load_template
=============

.. code-block:: python

    t = load_template('1F1T.pdb')
    t = load_template('1F1T.pdb','A')

load_alignment
==============

.. code-block:: python

    a = load_alignment('alignment_1F1T.fasta')
    
load_model
==========

.. code-block:: python

    m = load_model('1F1T.pdb')
    m = load_model('1F1T.pdb','A')
    m = load_model(biopy_struc, data_type='structure')
    m = load_model(biopy_struc[0].child_list[0], data_type='chain')


match_alignment_with_model
==========================

.. code-block:: python

    match_alignment_with_model(a,m)
    boolean = match_alignment_with_model(a,m)


match_template_with_alignment
=============================

.. code-block:: python

    match_template_with_alignment(t,a)
    boolean = match_template_with_alignment(t,a)


renumber_chain
==============

.. code-block:: python

    renumber_chain(m,'1')

remove_modification
===================

.. code-block:: python

    remove_modification(m['5'])
    remove_modification(t['5'], m)
    remove_modification(t['5'], m, '5A')


remove_all_modifications
========================

.. code-block:: python

    remove_all_modifications(m)


remove_mismatching_modifications
================================

.. code-block:: python

    remove_mismatching_modifications(t,a,m)


rotate_chi
==========

.. code-block:: python

    rotate_chi(m['5'], 90)


shrink_helix
============

.. code-block:: python

    shrink_helix(m, '2', '37', '5', '34')


write_logfile
=============

.. code-block:: python

    write_logfile()
    write_logfile('log.txt')


write_model
===========

.. code-block:: python

    write_model(m)
    write_model(m, 'output.pdb')
    write_model(m, 'output.pdb', 'log.txt')

write_fragment_candidates
=========================

.. code-block:: python

    candidates = find_fragment(m, '35', '38', 'CA', 20)
    write_fragment_candidates(candidates, 'my_candidates')

write_secstruc
==============

.. code-block:: python

    m = load_model('1F1T.pdb','A')
    write_secstruc(m, '1F1T_secondary_structure.vienna')


