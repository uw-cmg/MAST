###############################
Examples
###############################
***********************************************
Full example: defects, charges, NEB, phonons
***********************************************

Recipe::

    Recipe OptimizeWorkflow
    perfect_opt1 (lowmesh)
        perfect_opt2
            perfect_stat (static)
            {begin}
            inducedefect_<N> (inducedefect)
                defect_<N>_<Q>_opt1 (lowmesh_defect)
                    defect_<N>_<Q>_opt2 (defect_relax)
                        defect_<N>_<Q>_stat (static)
            {end}
    {begin}
    defect_<N>_<Q>_stat (static)
        phonon_<N>_<Q>_<P> (phonon)
            phonon_<N>_<Q>_<P>_parse (phononparse)
    {end}
    {begin}
    defect_<B>_<Q>_stat (static_to_neb), defect_<E>_<Q>_stat (static_to_neb)
        neb_<B-E>_<Q>_opt1 (neb_to_neb)
            neb_<B-E>_<Q>_opt2 (neb_to_nebstat)
                neb_<B-E>_<Q>_stat (nebstat_to_phonon)
        neb_<B-E>_<Q>_opt2 (neb_to_nebstat)
        neb_<B-E>_<Q>_stat (nebstat_to_phonon)
    {end}
    {begin}
    neb_<B-E>_<Q>_stat (nebstat_to_phonon)
        phonon_<B-E>_<Q>_<P> (phonon)
            phonon_<B-E>_<Q>_<P>_parse (phononparse)
    {end}

Input file::

    # Small demo for NEB workflow
    $mast
    system_name PhononNebTest
    $end

    $structure
    coord_type fractional

    begin elementmap
    X1 Al
    X2 Mg
    end

    begin lattice
    3.5 0 0
    0 3.5 0
    0 0 3.5
    end

    begin coordinates
    X1 0.0000000000 0.0000000000 0.0000000000
    X1 0.5000000000 0.5000000000 0.0000000000
    X1 0.0000000000 0.5000000000 0.5000000000
    X1 0.5000000000 0.0000000000 0.5000000000
    end

    $end

    $defects
    threshold 1e-4
    coord_type fractional

    begin int1
    interstitial 0.25 0.25 0.25 X2
    phonon host 0.0 0.5 0.5 0.5
    charge=-3,-2
    end

    begin int2
    interstitial 0.25 0.25 0.75 X2
    phonon host 0.0 0.0 0.0 0.5
    phonon int 0.25 0.25 0.75 0.5
    charge=-2,-2
    end

    begin int3
    interstitial 0.75 0.25 0.25 X2
    phonon host 0.0 0.0 0.0 0.5
    phonon int 0.75 0.25 0.25 0.5
    charge=-3,-3
    end

    $end

    $ingredients
    begin ingredients_global
    mast_nodes         1
    mast_multiplyencut 1.5
    mast_ppn           1
    mast_queue         default
    mast_exec          //share/apps/vasp5.2_cNEB
    mast_kpoints       2x2x2 M
    mast_xc            PBE
    isif 3
    ibrion 2
    nsw 191
    ismear 1
    sigma 0.2
    lwave False
    lcharg False
    prec Accurate
    mast_program   vasp
    mast_write_method            write_singlerun
    mast_ready_method            ready_singlerun
    mast_run_method              run_singlerun
    mast_complete_method         complete_singlerun
    mast_update_children_method  give_structure
    end

    begin inducedefect
    mast_write_method            no_setup
    mast_ready_method            ready_defect
    mast_run_method              run_defect
    mast_complete_method         complete_structure
    end

    begin lowmesh
    mast_kpoints 1x1x1 G
    end

    begin lowmesh_defect
    mast_kpoints 1x1x1 G
    isif 2
    end

    begin defect_relax
    isif 2
    end

    begin static
    ibrion -1
    nsw 0
    mast_multiplyencut 1.25
    mast_update_children_method give_structure
    end

    begin static_to_neb
    ibrion -1
    nsw 0
    mast_multiplyencut 1.25
    mast_update_children_method give_structure_and_energy_to_neb
    end

    begin phonon
    ibrion 5
    mast_write_method write_phonon_single
    mast_update_children_method give_phonon_single_forces_and_displacements
    end

    begin phononparse
    mast_program phon
    lfree .True.
    temperature 1173
    nd 3
    qa 11
    qb 11
    qc 11
    lsuper .False.
    mast_exec //home/tam/tammast/bin/phon_henry
    end

    begin neb_to_neb
    mast_kpoints 1x1x1 G
    mast_program   vasp_neb
    mast_write_method            write_neb
    mast_update_children_method  give_neb_structures_to_neb
    end

    begin neb_to_nebstat
    mast_program   vasp_neb
    mast_write_method            write_neb
    mast_update_children_method  give_neb_structures_to_neb
    end

    begin nebstat_to_phonon
    mast_program   vasp
    mast_write_method            write_neb_subfolders
    mast_ready_method            ready_neb_subfolders
    mast_run_method              run_neb_subfolders
    mast_complete_method         complete_neb_subfolders
    mast_update_children_method  give_saddle_structure
    end

    $end

    $neb
    begin int1-int2
    X2, 0.25 0.25 0.25, 0.25 0.25 0.75
    images 1
    phonon int 0.25 0.25 0.5 0.5
    phonon host 0.0 0.0 0.0 0.5
    end
    $end

    $recipe
    recipe_file phonon_test_neb.txt
    $end

*********************************************************
Small example: generic program (here, Genetic Algorithm)
*********************************************************
Recipe file::

    Recipe GenericTest
    generictest (generictest)

More lines could be added to the recipe, and more ingredient types (e.g. test1, test2, etc.), with minor modifications to the keywords given for each ingredient type.

Input file::

    $mast
    system_name GATest
    $end

    $structure
    #The structure actually does not make a difference for this
    #example, as it is not passed into any structure file.
    coord_type fractional
    begin lattice
    3.5 0 0
    0 3.5 0
    0 0 3.5
    end
    begin coordinates
    Al 0.0000000000 0.0000000000 0.0000000000
    end
    $end

    $ingredients
    begin ingredients_global
    mast_nodes         1
    mast_multiplyencut 1.5
    mast_ppn           1
    mast_queue         default
    mast_exec          //share/apps/vasp5.2_cNEB
    end

    begin generictest
    # need to add mastlib to python path to get lammps3.py
    # Amy's GAv14 is currently treated as closed-source
    mast_program                 None
    mast_exec                    python //home/tam/test_amy_GA/GAv14.py input.txt
    mast_complete_file           GAoutput.txt
    mast_complete_search         End of Execution
    mast_started_file            GAoutput.txt
    mast_copy_files //home/tam/tammast/test/gatest/SiC.tersoff //home/tam/tammast/test/gatest/cBulk.xyz
    mast_delimiter               =
    type  Defect
    atomlist  [('Si',0,28.0855,-5.3062),('C',4,12.011,-7.371)]
    filename   GAoutput
    nclust   5
    maxgen   5
    supercell   (3,3,3)
    SolidFile   cBulk.xyz
    SolidCell   [13.092,13.092,13.092]
    convergence_scheme   Max-Gen
    MUTPB  0.1
    mutation_options  ['Lattice_Alteration_small', 'Lattice_Alteration_Group', 'Rotation_geo']
    CALC_Method   LAMMPS
    pair_style  tersoff
    pot_file    SiC.tersoff
    LammpsMin   1e-25 1e-25 5000 10000
    keep_Lammps_files True
    Lmin_style   cg
    genealogy   True
    allenergyfile  True
    BestIndsList   True
    mast_write_method            write_singlerun
    mast_ready_method            ready_singlerun
    mast_run_method              run_singlerun
    mast_complete_method         complete_singlerun
    mast_update_children_method  give_structure
    end
    $end

    $recipe
    recipe_file generic_test.txt
    $end

