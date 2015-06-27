###################################
The Ingredients section
###################################

The ``$ingredients`` section contains a section for global ingredient keywords and then a section for each **ingredient type**.

Each ingredient type in the recipe should have a subsection denoted by ``begin <ingredient type>``.

Example ``$ingredients`` section::

    $ingredients

    begin ingredients_global
    keyword1 k1value1
    end

    begin ingredient_type1
    keyword1 k1value2
    keyword2 k2value1
    end

    begin ingredient_type2
    keyword2 k2value2
    end

    $end

Program-specific keywords such as VASP INCAR keywords are included in these sections. All other keywords are prefaced with ``mast_``. 

If there are no changes from the ingredients_global section, just add an empty subsection for that ingredient type::

    begin ingredient_type
    end

For a specific ingredient type, if a keyword is not specified in that ingredient type's subsection but is specified in the **ingredients_global** subsection, then, the value for that keyword will be taken from ingredients_global. 

*  In the example above, ``ingredient_type2`` would inherit ``keyword1 k1value1`` from ingredients_global.

========================================
Program-specific keywords
========================================

VASP keywords such as ``IBRION``, ``ISIF``, ``LCHARG``, ``LWAVE``, and so on, can be specified under each ingredient type in the ``$ingredients`` section of the input file.

Such program-specific keywords are only allowed if they are listed in the program-specific file located in the ``<MAST installation directory>/MAST/ingredients/programkeys/`` folder, for example, ``<MAST installation directory>/MAST/ingredients/programkeys/vasp_allowed_keywords.py``.

These program-specific keywords will be turned into uppercase keywords. The values will not change case, and should be given in the case required by the program. For example, ``lwave False`` will be translated into ``LWAVE False`` in the VASP INCAR file.

One exception for VASP keywords is the ``IMAGES`` keyword, which signals a nudged elastic band run, and should instead be set in the ``$neb`` section of the input file.

For VASP ingredients, please include ::

    lcharg False 
    lwave False 

in your ingredient global keywords in order to avoid writing the large VASP files CHGCAR and WAVECAR, unless you really need these files.

=================================
Special MAST keywords
=================================

Any keyword that starts with ``mast_`` is considered a special keyword utilized by MAST and will not be written into the VASP INCAR file or any custom input file.

-----------------------------------------
Submission script keywords
-----------------------------------------

The following queue submission keywords are platform-dependent and are used along to create the submission script (see :doc:`1_0_installation`).

**mast_exec**: The command used in the submission script to execute the program. Note that this is a specific command rather than the class of program, given in ``mast_program``, and it should include any MPI commands. ::

    mast_exec //opt/mpiexec/bin/mpiexec ~/bin/vasp_5.2

**mast_nodes**: The number of nodes requested.

**mast_ppn**: The number of processors per node requested.

**mast_queue**: The queue requested.

**mast_walltime**: The walltime requested, in whole number of hours

**mast_memory**: The memory per processor requested.

------------------------------------
MAST control flow keywords
------------------------------------

**mast_program**: Specify which program to run (``vasp``, ``vasp_neb``, or ``None`` for a generic program, are currently supported) ::

    mast_program vasp

*  This keyword must be in lowercase

**mast_frozen_seconds**: A number of seconds before a job is considered frozen, if its output file has not been updated within this amount of time. If not set, 21000 seconds is used.

**mast_auto_correct**: Specify whether mast should automatically correct errors.

*  The default is True, so if this keyword is set to True, or if this keyword is not specified at all, then MAST will attempt to find errors, automatically correct the errors, and resubmit the ingredient.
*  If set to False, MAST will attempt to find errors, then write them into a ``MAST_ERROR`` file in the recipe folder, logging both the error-containing ingredient and the nature of the error, but not taking any corrective actions. The recipe will be skipped in all subsequent MAST runs until the ``MAST_ERROR`` file is manually deleted by the user.


-----------------------------------
VASP-specific keywords
-----------------------------------

**mast_kpoints**: Specify k-point instructions in the form of kpoints along lattice vectors a, b, and c, and then a designation M for Monkhorst-Pack or G for Gamma-centered. :: 

    mast_kpoints = 3x3x3 G

*  Either this keyword or ``mast_kpoint_density`` is required for VASP calculations.

**mast_kpoint_density**: A number for the desired kpoint mesh density. 

*  Only works with ``mast_write_method`` of ``write_singlerun_automesh``
*  Either this keyword or ``mast_kpoints`` is required for VASP calculations.

**mast_pp_setup**: Specify which pseudopotential goes to which element::

    mast_pp_setup La=La Mn=Mn_pv O=O_s

**mast_xc**: Specify an exchange correlation functional; for VASP, follow the conventions of pymatgen (e.g. pw91, pbe)

*  This keyword is required for VASP calculations.

**mast_multiplyencut**: Specify a number with which to multiply the maximum ENCUT value of the pseudopotentials. Volume relaxations in VASP often take 1.5; otherwise 1.25 is sufficient.

*  Default is 1.5
*  If ``encut`` is given as a program keyword, then that value will be used and ``mast_multiplyencut`` should have no effect

**mast_setmagmom**: Specify a string to use for setting the initial magnetic moment. A short string will result in multipliers. For example, ``mast_setmagmom 1 5 1`` will produce ``2*1 2*5 8*1`` for a 12-atom unit cell with 2A, 2B, and 8C atoms. A string of the number of atoms in the ``POSCAR`` file will be printed as entered, for example, ``mast_setmagmom 1 -1 1 -1 1 -1 1 -1``.

**mast_charge**: Specify the charge on the system (total system)

*  -1 charge means the ADDITION of one electron. For example, O2- has two more electrons than O neutral. 
*  A positive charge is the REMOVAL of electrons. For example, Na+ with a +1 charge has one FEWER electron than Na neutral.

**mast_coordinates**: For a non-NEB calculation, allows you to specify a single POSCAR-type of CIF structure file which corresponds to the relaxed fractional coordinates at which you would like to start this ingredient. ONLY the coordinates are used. The lattice parameters and elements are given by the $structure section of the input file. The coordinates must be fractional coordinates. ::

    mast_coordinates POSCAR_initialize

*  For an NEB calculation, use a comma-delimited list of poscar files corresponding to the correct number of images. Put no spaces between the file names. Example for an NEB with 3 intermediate images::
    
    mast_coordinates POSCAR_im1,POSCAR_im2,POSCAR_im3

*  The structure files must be found in the directory from which the input file is being submitted when initially inputting the input file (e.g. the directory you are in when you run ``mast -i test.inp``); once the ``input.inp`` file is created in the recipe directory, it will store a full path back to these poscar-type files.

*  This keyword cannot be used with programs other than VASP, cartesian coordinates, and special ingredients like inducedefect-type ingredients, whose write or run methods are different.

---------------------------------
Structure manipulation keywords
---------------------------------
**mast_strain**: Specify three numbers for multiplying the lattice parameters a, b, and c. Only works with ``mast_run_method`` of ``run_strain`` ::

    mast_strain 1.01 1.03 0.98 

This example will stretch the lattice along lattice vector a by 1%, stretch the lattice along lattice vector b by 3%, and compress the lattice along lattice vector c by 2%

---------------------------
mast_xxx_method keywords
---------------------------

The following keywords have individual sections:

**mast_write_method**: Specifies what the ingredient should write out before running (e.g., create the INCAR)

**mast_ready_method**: Specifies how MAST can tell if the ingredient is ready to run (often, in addition to writing its own files, an ingredient must also wait for data from its parent ingredient(s)). 

**mast_run_method**: Specifies what MAST should do to run the ingredient (e.g. submit a submission script to a queue, or perform some other action)

**mast_complete_method**: Specifies how MAST can tell if the ingredient is considered complete

**mast_update_children_method**: Specifies what information an ingredient passes on to its children, and how it does so.

   
Specific available values for each keyword are given in the accompanying sections, and require no arguments, e.g.::

    mast_write_method write_singlerun

They depend on having an appropriate program set in ``mast_program``.

However, you may choose to specify arguments where available, e.g.::
    
    mast_complete_method file_has_string myoutput "End of Execution"

You may also choose to specify multiple methods. These methods will be performed in the order listed. For ``mast_ready_method`` or ``mast_complete_method``, all methods listed must return True in order for the ingredient to be considered ready or complete, respectively. 
Use a semicolon to separate out the methods::

    mast_complete_method file_has_string myoutput "End of Execution"; file_exists Parsed_Structures

In the example above, the file "myoutput" must exist and contain the phrase "End of Execution", and the file "Parsed_Structures" must exist, in order for the ingredient to be considered complete.

Update-children methods will always get the child name appended as the end of the argument string. For example, ::

    mast_update_children_method copy_file EndStructure BeginStructure

will copy the file EndStructure of the parent ingredient folder to a new file BeginStructure in the child ingredient folder. There is no separate argument denoting the child ingredient.

All arguments are passed as strings. Arguments in quotation marks are kept together.

Some common open-ended methods are:

*  **file_exists <filename>**

*  **file_has_string <filename> <string>**

*  **copy_file <filename> <copy_to_filename>**

*  **softlink_file <filename> <softlink_to_filename>**

*  **copy_fullpath_file <full path file name> <copy_to_filename>**: This method is for copying some system file like //home/user/some_template, not an ingredient-specific file

*  **write_ingred_input_file <filename> <allowed file> <uppercase keywords> <delimiter>**: The allowed file specifies an allowed keywords file name in ``<MAST installation directory>/MAST/ingredients/programkeys``. 

    *  Use "all" to put any non-mast keywords into the input file. 
    *  Use 1 to uppercase all keywords, or 0 to leave them as entered. 
    *  Leave off the delimiter argument in order to use a single space. 
    *  Examples::
    
        write_ingred_input_file input.txt all 0 =
        write_ingred_input_file input.txt phon_allowed_keys.py 1

*  **no_setup**: Does nothing. Useful when you want to specifically specify doing nothing.

*  **no_update**: Does nothing (but, does accept the child name it is given). Useful when you want to specify doing nothing for a child update step.

*  **run_command: <command string, including all arguments>**: This method allows you to run a python script. 

    *  The python script may take in only string-based arguments
    *  Please stick to common text characters. 
    *  Example:: 
    
        mast_run_method run_command "//home/user/myscripts/my_custom_parsing.py 25 0.01"

    *  In the example above, the numbers 25 and 0.01 will actually be passed into sys.argv as a string. 
    *  This method is intended to allow you to run short custom scripts of your own creation, particularly for ``mast_write_method`` when setting up your ingredient.
    *  For long or complex execution steps where you want the output tracked separately, do not use this method. Instead, do the following in order to get your script submitted to the queue:
        #  Use ``write_submit_script`` in your ``mast_write_method``, along with any other write methods
        #  Use ``mast_run_method run_singlerun``
        #  Put your script in the ``mast_exec`` keyword

-----------------------------------------
mast_write_method keyword values
-----------------------------------------

**write_singlerun**

*  Write files for a single generic run.
*  Programs supported: vasp 

**write_singlerun_automesh**

*  Write files for a single generic run.
*  Programs supported: vasp
*  Requires the ``mast_kpoint_density`` ingredient keyword

**write_neb**

*  Write an NEB ingredient. This method writes interpolated images to the appropriate folders, creating 00/01/.../0N directories and uses linear interpolation between images.
*  Programs supported: vasp

**write_pathfinder_neb**

*  Write an NEB ingredient using a charge-density-based pathfinding method.
*  Programs supported: vasp
*  This method takes the argument of an ingredient name for the ingredient from which to take the charge density.

    * The ingredient name must be fully specified, e.g. no <S>, <N>, etc. tags. 
    * The ingredient must have a CHGCAR file written. A gamma-point calculation is sufficient for this purpose.
    * The ingredient should have both endpoints removed. For example, for vacancy migration, the ingredient should have neither a vacancy at the initial position, nor a vacancy at the final position. However, it should have all other non-migrating defects that are common to both the initial and final state.

*  This method only works with ``use_structure_index True`` in the ``$structure`` section of the input file.

**write_neb_subfolders**

*  Write static runs for an NEB, starting from a previous NEB, into image subfolders 01 to 0(N-1).
*  Programs supported: vasp

**write_phonon_single**

*  Write files for a phonon run.
*  Programs supported: vasp

**write_phonon_multiple**

*  Write a phonon run, where the frequency calculation for each atom and each direction is a separate run, using selective dynamics. CHGCAR and WAVECAR must have been given to the ingredient previously; these files will be softlinked into each subfolder.
*  Programs supported: vasp

-----------------------------------------
mast_ready_method keyword values
-----------------------------------------

**ready_singlerun**

*  Checks that a single run is ready to run
*  Programs supported: vasp (either NEB or regular VASP run), phon

**ready_defect**

*  Checks that the ingredient has a structure file
*  Programs supported: vasp

**ready_neb_subfolders**

*  Checks that each 01/.../0(N-1) subfolder is ready to run as its own separate calculation, following the ready_singlerun criteria for each folder
*  This method is used for NEB static calculations rather than NEB calculations themselves.

**ready_subfolders**
*  Checks that each subfolder is ready to run, following the ready_singlerun criteria.
*  Generic
*  This method is used for calculations whose write method includes subfolders, and where each subfolder is a calculation, as in ``write_phonon_multiple``.

----------------------------------
mast_run_method keyword values
----------------------------------

**run_defect**

*  Create a defect in the structure; not submitted to queue
*  Generic
*  Requires the ``$defects`` section in the input file (see :doc:`3_1_5_defects`).

**run_singlerun**

*  Submit a run to the queue.
*  Generic

**run_neb_subfolders**

*  Run each 01/.../0(N-1) subfolder as run_singlerun
*  Generic

**run_subfolders**

*  Run each subfolder as run_singlerun
*  Generic

**run_strain**

*  Strain the structure; not submitted to queue
*  Generic
*  Requires the ``mast_strain`` ingredient keyword

**run_scale**

*  Scale the structure (e.g. a 2-atom unit cell scaled by 2 becomes a 16-atom supercell)
*  Requires the ``$scaling`` subsection in the input file (see :doc:`3_1_1_structure`).
*  Must not be run on the starting ingredient.


------------------------------------
mast_complete_method keyword values
------------------------------------

**complete_singlerun**

*  Check if run is complete
*  Programs supported: vasp 
*  Note that for VASP, the phrase ``reached required accuracy`` is checked for, as well as a ``User time`` in seconds. The exceptions are:

    *  NSW of 0, NSW of -1, or NSW not specified in the ingredients section keywords is taken as a static calculation, and .EDIFF is reached. is checked instead of .reached required accuracy.
    *  IBRION of -1 is taken as a static calculation, and .EDIFF is  reached. is checked instead of .reached required accuracy.
    *  IBRION of 0 is taken as an MD calculation, and only user time is checked
    *  IBRION of 5, 6, 7, or 8 is taken as a phonon calculation, and only user time is checked

**complete_neb_subfolders**

*  Check if all NEB subfolders 01/.../0(N-1) are complete, according to complete_singlerun criteria.
*  This method is not for checking the completion of NEBs! An NEB ingredient should have ``mast_program vasp_neb`` and ``mast_complete_method complete_singlerun``.
*  An NEB static calculation, or a static calculation for each image, would use this keyword as ``mast_complete_method complete_neb_subfolders`` but have ``mast_program vasp`` instead of vasp_neb.

**complete_subfolders**

*  Check if all subfolders are complete, according to complete_singlerun criteria.
*  Generic

**complete_structure**

*  Check if run has an output structure file written
*  Programs supported: vasp (looks for CONTCAR)

--------------------------------------------
mast_update_children_method keyword values
--------------------------------------------

**give_structure**

*  Forward the relaxed structure
*  Programs supported: vasp (CONTCAR to POSCAR)

**give_structure_and_energy_to_neb**

*  Forward the relaxed structure and energy files
*  Programs supported: vasp (CONTCAR to POSCAR, and copy over OSZICAR)

**give_neb_structures_to_neb**

*  Give NEB output images structures as the starting point image input structures in another NEB
*  Programs supported: vasp (01/.../0(N-1) CONTCAR files will be the child NEB ingredient.s starting 01/.../0(N-1) POSCAR files.

**give_saddle_structure**

*  Forward the highest-energy structure of all subfolder structures
*  Programs supported: vasp

-------------------------------
Example Ingredients section
-------------------------------

Here is an example ingredients section::

    $ingredients
    begin ingredients_global
    mast_program    vasp
    mast_nodes      1
    mast_multiplyencut 1.5
    mast_ppn        1
    mast_queue      default
    mast_exec       mpiexec //home/mayeshiba/bin/vasp.5.3.3_vtst_static
    mast_kpoints    2x2x2 M
    mast_xc PW91
    isif 2
    ibrion 2
    nsw 191
    ismear 1
    sigma 0.2
    lwave False
    lcharg False
    prec Accurate
    mast_program   vasp
    mast_write_method           write_singlerun
    mast_ready_method           ready_singlerun
    mast_run_method             run_singlerun
    mast_complete_method        complete_singlerun
    mast_update_children_method  give_structure
    end

    begin volrelax_to_singlerun
    isif 3
    end

    begin singlerun_to_phonon
    ibrion -1
    nsw 0
    mast_update_children_method  give_structure_and_restart_files
    mast_multiplyencut 1.25
    lwave True
    lcharge True
    end

    begin inducedefect
    mast_write_method           no_setup
    mast_ready_method           ready_defect
    mast_run_method             run_defect
    mast_complete_method        complete_structure
    end

    begin singlerun_vac1
    mast_coordinates            POSCAR_vac1
    end

    begin singlerun_vac2
    mast_coordinates            POSCAR_vac2
    end

    begin singlerun_to_neb
    ibrion -1
    nsw 0
    mast_update_children_method  give_structure_and_energy_to_neb
    lwave True
    lcharge True
    end

    begin neb_to_neb_vac1-vac2
    mast_coordinates            POSCAR_nebim1,POSCAR_nebim2,POSCAR_nebim3
    mast_write_method           write_neb
    mast_update_children_method  give_neb_structures_to_neb
    mast_nodes                  3
    mast_kpoints                1x1x1 G
    ibrion 1
    potim 0.5
    images 3
    lclimb True
    spring -5
    end

    begin neb_to_neb_vac1-vac3
    mast_coordinates            POSCAR_nebim1_set2,POSCAR_nebim2_set2,POSCAR_nebim3_set2
    mast_write_method           write_neb
    mast_update_children_method  give_neb_structures_to_neb
    mast_nodes                  3
    mast_kpoints                1x1x1 G
    ibrion 1
    potim 0.5
    images 3
    lclimb True
    spring -5
    end

    begin neb_to_nebstat
    mast_write_method           write_neb
    mast_update_children_method  give_neb_structures_to_neb
    mast_nodes                  3
    ibrion 1
    potim 0.5
    images 3
    lclimb True
    spring -5
    end

    begin nebstat_to_nebphonon
    ibrion -1
    nsw 0
    mast_write_method           write_neb_subfolders
    mast_ready_method           ready_neb_subfolders
    mast_run_method             run_neb_subfolders
    mast_complete_method        complete_neb_subfolders
    mast_update_children_method  give_saddle_structure
    end

    begin phonon_to_phononparse
    mast_write_method           write_phonon_multiple
    mast_ready_method           ready_subfolders
    mast_run_method             run_subfolders
    mast_complete_method        complete_subfolders
    mast_update_children_method  give_phonon_multiple_forces_and_displacements
    ibrion 5
    nfree 2
    potim 0.01
    istart 1
    icharg 1
    end
    
    $end

