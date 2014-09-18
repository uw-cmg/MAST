###################################
Input File
###################################

**********************************
Introduction to the Input File
**********************************

The MAST program is driven by two main files: an input file which contains all the various keywords required for setting up the recipe (i.e. workflow), and a :doc:`Recipe Template <4_0_recipe>` which organizes all the ingredients (i.e. calculations) in the recipe. In this section, we will discuss the input file

The input file contains several sections and subsections.
Bounds of sections are denoted by ``$sectionname`` and ``$end``.
Bounds of subsections within a section are denoted by ``begin subsectionname`` and ``end``.

*Comments in the input file are allowed only as separate lines starting with #. A comment may not be appended to a line.*

Example of the ``$structure`` section, with three subsections, **elementmap**, **coordinates**, and **lattice**::

    $structure
    coord_type fractional
    
    begin elementmap
    X1 Ga
    X2 As
    end
    
    begin coordinates
    X1 0.000000 0.000000 0.000000
    X1 0.500000 0.500000 0.000000
    X2 0.250000 0.250000 0.250000
    X2 0.750000 0.750000 0.250000
    end
    
    begin lattice
    6.0 0.0 0.0
    0.0 6.0 0.0
    0.0 0.0 6.0
    end
    $end

Each section is described in detail below.

************************
The MAST section
************************
The ``$mast`` section contains this keyword:

*  system_name: Specify a single descriptive word here, like EpitaxialStrain. This keyword will become part of the recipe directory.s name and allow you to spot the recipe in the ``$MAST_SCRATCH`` directory::

    system_name EpitaxialStrain


*****************************
The Structure section
*****************************

The ``$structure`` section contains the coordinate type, coordinates, and lattice, or, optionally, the name of a structure file (either CIF or VASP POSCAR-type).

====================================
Structure by file
====================================

Using the keyword ``posfile``, a VASP POSCAR-type file or a CIF file can be inserted here in this section::

    $structure
    posfile POSCAR_fcc
    $end

The file should be located in the same directory as the input file.

A CIF file should end with .cif.

A POSCAR-type filename must start with ``POSCAR_`` or ``CONTCAR_`` in order for pymatgen to recognize it. The elements will be obtained from the POSCAR unless you also have a POTCAR in the directory, in which case, check your output carefully because the elements might be given by the POTCAR instead, no matter what elements are written in the POSCAR file.

====================================
Structure by specification
====================================

To specify a structure, use the following subsections:

**coord_type**: This keyword specifies fractional or cartesian coordinates. Only fractional coordinates have been thoroughly tested with most MAST features.

**lattice**: The lattice subsection specifies lattice basis vectors on a cartesian coordinate system.

**elementmap**: The elementmap subsection allows you to create a generic lattice and interchange other elements onto it. This is useful when looping over other elements (discussed later).

The elementmap subsection works in conjunction with the coordinates subsection.

**coordinates**: The coordinates subsection specifies the coordinates in order. 

Fractional coordinates are fractional along each lattice basis vector, e.g. .0.5 0 0. describes a position 0.5 (halfway) along the first lattice basis vector.

Each fractional coordinate must be preceded by either an element symbol or an X# symbol corresponding to the symbols assigned in the elementmap section.


Example::
    
    begin $structure

    coord_type fractional    

    begin lattice
    6.0 0.0 0.0
    0.0 6.0 0.0
    0.0 0.0 6.0
    end

    begin elementmap
    X1 Ga
    X2 As
    end
    
    begin coordinates
    X1 0.000000 0.000000 0.000000
    X1 0.500000 0.500000 0.000000
    X1 0.000000 0.500000 0.500000
    X1 0.500000 0.000000 0.500000
    X2 0.250000 0.250000 0.250000
    X2 0.750000 0.750000 0.250000
    X2 0.250000 0.750000 0.750000
    X2 0.750000 0.250000 0.750000
    end
    
    $end

*************************
The Ingredients section
*************************

The ``$ingredients`` section contains a section for global ingredient keywords and then a section for each ingredient type. 

Program-specific keywords such as VASP INCAR keywords are included in these sections. All other keywords are prefaced with ``mast_``. 

Each ingredient type in the recipe should have a subsection denoted by ::

    begin ingredient_type
    (keywords here)
    end

even if there are no keywords within that section, in which case the ``end`` line directly follows the ``begin`` line.

========================================
Ingredients that are VASP calculations
========================================

VASP keywords such as ``IBRION``, ``ISIF``, ``LCHARG``, ``LWAVE``, and so on, can be specified under each ingredient type in the ``$ingredients`` section of the input file.

Such program-specific keywords are only allowed if they are listed in the program-specific file located in the ``$MAST_INSTALL_PATH/MAST/ingredients/programkeys/`` folder, for example, ``$MAST_INSTALL_PATH/MAST/ingredients/programkeys/vasp_allowed_keywords.py``.

These program-specific keywords will be turned into uppercase keywords. The values will not change case, and should be given in the case required by the program. For example, ``lwave False`` will be translated into ``LWAVE False`` in the VASP INCAR file.

One exception for VASP keywords is the ``IMAGES`` keyword, which signals a nudged elastic band run, and should instead be set in the ``$neb`` section of the input file.

For VASP ingredients, please include ::

    lcharg False 
    lwave False 

in your ingredient global keywords in order to avoid writing the large VASP files CHGCAR and WAVECAR, unless you really need these files.


Any keyword that starts with ``mast_`` is considered a special keyword utilized by MAST and will not be written into the VASP INCAR file.

===================================
Special MAST ingredient keywords:
===================================

Some of these special MAST keywords are only appropriate for VASP calculations.

**mast_program**: Specify which program to run (``vasp``, ``vasp_neb``, ``phon``, or ``None`` for a generic program, are currently supported) ::

    mast_program vasp

*  This keyword must be in lowercase

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

**mast_strain**: Specify three numbers for multiplying the lattice parameters a, b, and c. Only works with ``mast_run_method`` of ``run_strain`` ::

    mast_strain 1.01 1.03 0.98 

This example will stretch the lattice along lattice vector a by 1%, stretch the lattice along lattice vector b by 3%, and compress the lattice along lattice vector c by 2%

**mast_scale**: A number for which to scale all dimensions of a supercell. Only works with ``mast_run_method`` of ``run_scale`` or ``run_scale_defect``

**mast_frozen_seconds**: A number of seconds before a job is considered frozen, if its output file has not been updated within this amount of time. If not set, 21000 seconds is used.

**mast_auto_correct**: Specify whether mast should automatically correct errors.

*  The default is True, so if this keyword is set to True, or if this keyword is not specified at all, then MAST will attempt to find errors, automatically correct the errors, and resubmit the ingredient.
*  If set to False, MAST will attempt to find errors, then write them into a ``MAST_ERROR`` file in the recipe folder, logging both the error-containing ingredient and the nature of the error, but not taking any corrective actions. The recipe will be skipped in all subsequent MAST runs until the ``MAST_ERROR`` file is manually deleted by the user.

The following keyword is used only for generic programs (not VASP, PHON, or any other named programs). 

**mast_started_file**: A file name in the ingredient directory whose presence signals that the ingredient run has been started. ::

    mast_started_file        GAoutput.txt

The following queue-submission keywords are platform dependent and are used along to create the submission script:

**mast_exec**: The command used in the submission script to execute the program. Note that this is a specific command rather than the .class. of program, given in ``mast_program``, and it should include any MPI commands. ::

    mast_exec //opt/mpiexec/bin/mpiexec ~/bin/vasp_5.2

**mast_nodes**: The number of nodes requested.

**mast_ppn**: The number of processors per node requested.

**mast_queue**: The queue requested.

**mast_walltime**: The walltime requested, in whole number of hours

**mast_memory**: The memory per processor requested.


The following keywords have individual sections:

**mast_write_method**: The .write. method, which specifies files the ingredient should write out before running (e.g., create the INCAR) 

**mast_ready_method**: The .ready. method, which specifies how MAST can tell if the ingredient is ready to run (often, in addition to writing its own files, an ingredient must also wait for data from its parent ingredient(s)). 

**mast_run_method**: The .run. method, which specifies what MAST should do to run the ingredient (e.g. submit a submission script to a queue, or perform some other action)

**mast_complete_method**: The .complete. method, which specifies how MAST can tell if the ingredient is considered complete

**mast_update_children_method**: the .update children. method, which specifies what information an ingredient passes on to its children, and how it does so.

.. _important_notes:
   
--------------------------------------------------
Important notes on using mast_xxx_method keywords
--------------------------------------------------
Specific available values for each keyword are given in the accompanying sections, and require no arguments, e.g.::

    mast_write_method write_singlerun

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

*  **write_ingred_input_file <filename> <allowed file> <uppercase keywords> <delimiter>**: The allowed file specifies an allowed keywords file name in ``$MAST_INSTALL_PATH/MAST/ingredients/programkeys``. 

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
        *  Some useful scripts are found in ``$MAST_INSTALL_PATH/tools`` and described in :ref:`6_0_tools`

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

*  Write an NEB ingredient. This method writes interpolated images to the appropriate folders, creating 00/01/.../0N directories.
*  Programs supported: vasp

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
*  Requires the ``$defects`` section in the input file.

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
*  Generic
*  Requires the ``mast_scale`` ingredient keyword, and must not be run on the starting ingredient (for VASP, the ingredient must already have been given a smaller POSCAR file, like the POSCAR for a 2-atom unit cell)

**run_scale_defect**

*  Scale the structure and defect it (e.g. a single defect at 0.5 0.5 0.5 in the original structure becomes a single defect at 0.25 0.25 0.25 in the structure scaled by 2)
*  Generic
*  Requires the ``mast_scale`` ingredient keyword, and must not be run on the starting ingredient

----------------------------------
mast_complete_method keyword values
----------------------------------

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

The following keywords are deprecated. Please use the generic methods in :ref:`important_notes` instead.

give_structure_and_restart_files (same as give_structure_and_restart_files_softlinks)

*  Forward the relaxed structure and additional files
*  Programs supported: vasp (CONTCAR to POSCAR, and softlinks to parent.s WAVECAR and CHGCAR files)

give_structure_and_restart_files_full_copies

*  Forward the relaxed structure and additional files
*  Programs supported: vasp (CONTCAR to POSCAR, and full copies of parent.s WAVECAR and CHGCAR files)

give_structure_and_charge_density_full_copy

*  Forward the relaxed structure and charge density file; copies the file
*  Programs supported: vasp (CONTCAR to POSCAR, and copy over CHGCAR)

give_structure_and_charge_density_softlink

*  Forward the relaxed structure and charge density file as a softlink
*  Programs supported: vasp (CONTCAR to POSCAR, and softlink to CHGCAR)

give_structure_and_wavefunction_full_copy
*  Forward the relaxed structure and wavefunction file; copies the file
*  Programs supported: vasp (CONTCAR to POSCAR, and copy over WAVECAR)

give_structure_and_wavefunction_softlink

*  Forward the relaxed structure and wavefunction file as a softlink
*  Programs supported: vasp (CONTCAR to POSCAR, and softlink to WAVECAR)

----------------------------------
Custom mast_xxx_method keywords
----------------------------------
You may also choose to write your own methods, in addition to any of the methods above.

Place these methods in a file in the directory ``$MAST_INSTALL_PATH/customlib``, structured like the file ``$MAST_INSTALL_PATH/customlib/customchopingredient.py``

*  Please inherit from either ChopIngredient or BaseIngredient.
*  Name the method(s) something unique (e.g. not found in either ChopIngredient or BaseIngredient)
*  You will have access to the ingredient directory name at ``self.keywords['name']`` as well as ingredient keywords at ``self.keywords['program_keys']``.
*  The method may also take in up to 3 string-based arguments.
*  In the input file, designate your custom method as classname.methodname followed by any arguments, for example, ``mast_write_method MyChopClass.write_complex_file superfile``


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
    
    begin phononparse
    mast_program                phon
    lfree .True.
    temperature 273
    ptemp 10 110
    nd 3
    qa 11
    qb 11
    qc 11
    lnosym .True.
    ldrift .False.
    lsuper .False.
    mast_exec $MAST_INSTALL_PATH/bin/phon_henry
    mast_multiplyencut 1.25
    end
    
    $end


********************
The Recipe section
********************

The ``$recipe`` section contains the recipe template to be used. More information on the recipe template is given in :doc:`4_0_recipe`.::

    $recipe
    perfect_opt1 (relax)
        perfect_stat (static)
        {begin}
        inducedefect_<N> (inducedefect)
            defect_<N>_opt1 (defect_relax)
                defect_<N>_stat (defect_static)
        {end}
    $end

*******************************
The Defects section (optional)
*******************************

The ``$defects`` section includes the defect type of vacancy, interstitial, substitution, or antisite (which is the same as substitution), the defect coordinates, and the defect element symbol.

*  Note that if an ``elementmap`` subsection is given in the ``$structure`` section of the input file, the mapped designations ``X1``, ``X2``, and so on can be given instead of an element symbol.

The ``coord_type`` keyword specifies fractional or cartesian coordinates for the defects.

The ``threshold`` keyword specifies the absolute threshold for finding the defect coordinate, since relaxation of the perfect structure may result in changed coordinates.

Example ``$defects`` section::

    $defects

    coord_type fractional
    threshold 1e-4

    vacancy 0 0 0 Mg
    vacancy 0.5 0.5 0.5 Mg
    interstitial 0.25 0.25 0 Mg
    interstitial 0.25 0.75 0 Mg
    
    $end

The above section specifies 4 point defects (2 vacancies and 2 interstitials) to be applied separately and independently to the structure. When combined with the correct :doc:`recipe <4_0_recipe>`, four separate ingredients, each containing one of the defects above, will be created.

Multiple point defects can be also grouped together as a combined defect within a .begin/end,. with a label after the .begin,. such as::

    $defects
    
    coord_type fractional
    threshold 1e-4
    
    begin doublevac
    vacancy 0.0 0.0 0.0 Mg
    vacancy 0.5 0.5 0.5 Mg
    end
    
    interstitial 0.25 0.25 0 Mg
    interstitial 0.25 0.75 0 Mg
    
    $end

In this case, there will be three separate .defect. ingredients: one ingredient with two vacancies together (where the defect group is labeled .doublevac.), one interstitial, and another interstitial.

Charges can be specified as ``charge=0,10``, where a comma denotes the lower and upper ranges for the charges.

Let's say we want a Mg vacancy with charges from 0 to 3 (0, 1, 2, and 3)::

    vacancy 0 0 0 Mg charge=0,3

Let.s say we want a dual Mg vacancy with a charge from 0 to 3 and labeled as Vac@Mg-Vac@Mg::

    begin Vac@Mg-Vac@Mg
    vacancy 0.0 0.0 0.0 Mg
    vacancy 0.5 0.5 0.5 Mg
    charge=0,3
    end

For a single defect, charges and labels can be given at the same time:

Let's say we have a Mg vacancy with charges between 0 and 3, and we wish to label it as Vac@Mg::

    vacancy 0.0 0.0 0.0 Mg charge=0,3 label=Vac@Mg

The charge and label keywords are interchangeable, i.e. we could also have typed::

    vacancy 0 0 0 Mg label=Vac@Mg charge=0,3

If you use charges in the defects section like this, then you should use a :doc:`recipe <4_0_recipe>` template with a free-form defect_<N>_<Q> format. 

=====================
Phonons for defects
=====================

Phonon calculations are described by a *phonon center site* coordinate and a *phonon center radius* in Angstroms. Atoms within the sphere specified by these two values will be included in phonon calculations.

For VASP, this inclusion takes the form of selective dynamics T T T for the atoms within the sphere, and F F F otherwise, in a phonon calculation (IBRION = 5, 6, 7, 8)

If the phonon center radius is 0, only the atom found at the phonon center site point will be considered.

To use phonons in the defects section, use the subsection keyword .phonon. followed by a label for the phonon, the fractional coordinates for the phonon center site, a float value for the phonon center radius, and an optional float value for the tolerance-matching threshold for matching the phonon center site (if this last value is not specified, 0.1 is used). Multiple separate phonon calculations may be obtained for each defect, for example::

    begin int1
    interstitial 0.25 0.25 0.25 X2
    phonon host3 0.3 0.3 0.4 2.5 0.01
    phonon solute 0.1 0.1 0.2 0.5
    end

In the example above, *host3* is the label for the phonon calculation where (0.3, 0.3, 0.4) is the coordinate for the phonon center site, and 2.5 Angstroms is the radius for the sphere inside which to consider atoms for the phonon calculation. Points within 0.01 of fractional coordinates will be considered for matching the phonon center site. 

In the example above, *solute* is the label for the phonon calculation bounded within a 0.5 Angstrom radius centered at (0.1, 0.1, 0.2) in fractional coordinates. As no threshold value was given, points within 0.1 (default) of fractional coordinates will be considered for matching the phonon center site.

The recipe template file for phonons may include either the explicit phonon labels and their charge and defect label, or <N>_<Q>_<P> (defect label _ charge label _ phonon label).

Because phonons are cycled with the defects, a new parent loop must be provided for the phonons, for example::

    {begin}
    defect_<N>_<Q>_stat (static)
        phonon_<N>_<Q>_<P> (phonon)
            phonon_<N>_<Q>_<P>_parse (phononparse)
    {end}

*********************************
The chemical potentials section
*********************************

The $chemical_potentials section lists chemical potentials, used for defect formation energy calculations using the defect formation energy tool.
Currently, chemical potentials must be set ahead of time. Each chemical potential set may be labeled. ::

    $chemical_potentials
    
    begin Ga rich
    Ga -3.6080
    As -6.0383
    Bi -4.5650
    end
    
    begin As rich
    Ga -4.2543
    As -5.3920
    Bi -4.5650
    end
    
    $end

********************
The NEB section
********************

The ``$neb`` section includes a list of nudged-elastic-band hops. Each neb hop should be a subsection labeled with the starting and ending .defect group. as specified in the ``$defects`` section, and then also indicate the movement of elements, and their closest starting and ending positions. These explicit positions disambiguate between possible interpolations.

*  Note that if an ``elementmap`` subsection is given in the ``$structure`` section of the input file, the mapped designations ``X1``, ``X2``, and so on can be given instead of an element symbol.

Again, the ``$neb`` section is tied to specific defect labels. The NEB ingredients must be able to find defects or defect groups with those labels.

The ``images`` keyword specifies the number of intermediate images, which must currently be the same in all NEBs in the recipe. 

Phonons may be specified within each NEB grouping, as in the defects section. The presumed saddle point in an NEB is usually taken; use the ``mast_update_children give_saddle_structure`` to give that saddle point structure to the phonon calculation. If, in an NEB, the frequencies for the moving atom are desired for the phonon calculations, and if that atom is anticipated to pass from fractional coordinate 0 0 0 to fractional coordinate 0.5 0 0, then the phonon_center_site should be 0.25 0 0 (assuming a straight path), and the phonon_center_radius is probably about 1 Angstrom. 

Example defect and NEB section together::

    $defects
    
    coord_type fractional
    threshold 1e-4
    
    vacancy 0.0 0.0 0.0 Mg label=vac1
    vacancy 0.0 0.5 0.5 Mg label=vac2
    interstitial 0.25 0.0 0.0 Al label=int1
    interstitial 0.0 0.25 0.0 Al label=int2
    
    $end
    
    $neb
    
    begin vac1-vac2
    images 1
    Mg, 0 0 0, 0 .5 0.5
    end
    
    begin int1-int2
    Al, 0.25 0 0, 0 0.25 0
    images 3
    phonon movingatom 0.125 0.125 0.0 1.0
    end
    
    $end

**************************************************************
Creating several input files at once: the looped input file
**************************************************************

One input file may be able to spawn several nearly-identical input files, which differ in small ways.

======================
Independent loops
======================

The special looping keyword ``indeploop`` may be used to signify a line which indicates that spawned input files should cycle through these values. ::

   indeploop mast_xc (pw91, pbe)

In this example, two input files will be created. One input file will contain the line ``mast_xc pw91``. The other input file will contain the line ``mast_xc pbe``.

*  Any text within parentheses and separated by a comma will be looped. 

*  Lines which normally include commas, like the ``charge`` line in the ``$defects`` section, or the ``mast_coordinates`` keyword for an NEB, may not be looped.

*  This keyword may only be used once on a line.

If there is more than one ``indeploop`` keyword in the input file, a combinatorial spawn of input files will be created.

For example, this excerpt would generate four input files: one with iron using pw91, one with iron using pbe, one with copper using pw91, and one with copper using pbe::
    
    $structure
    begin elementmap
    indeploop X1 (Fe, Cu)
    end
    ...
    $end
        
    $ingredients
    begin ingredients_global
    indeploop mast_xc (pw91, pbe)
    ...
    end
    $end

=============================
Dependent, or pegged, loops
=============================

Sometimes looped lines should really be looped together at the same time, rather than with each value looped over each other value.
 
For example, if you want to create a single input file, but signify that it should be copied into three input files, one for each element, but with different GGA+U U-values, you would use a pegged loop like this::

    $structure
    begin elementmap
    pegloop1 X1 (Es, Fm, Md)
    end
    ...
    $end
    
    $ingredients
    begin ingredients_global
    pegloop1 ldauu (5.3, 6.5, 8.0)
    ...
    end
    $end

In this case, three input files will be created. In the first input file, Es will be paired with a U-value of 5.3. In the second input file, Fm will be paired with a U-value of 6.5. In the third input file, Md will be paired with a U-value of 8.0.

There are two pegged loops allowed, specified by ``pegloop1`` and ``pegloop2``. 

Each pegged loop and independent loop will be combinatorially combined. For example, if a separate line ``indeploop mast_xc (pw91, pbe)`` were included in the ``ingredients_global`` subsection above, then six input files would be created: one pw91 and one pbe input file for Es with +U 5.3, another pair for Fm, and another pair for Mn.

In the example below, four input files would be created, corresponding to four different lattices:
*  [(6.0,0.0,0.0),(0.0,6.0,0.0),(0.0,0.0,2.0)]
*  [(6.0,0.0,0.0),(0.0,6.0,0.0),(0.0,0.0,3.0)]
*  [(4.0,0.0,0.0),(0.0,4.0,0.0),(0.0,0.0,2.0)], and 
*  [(4.0,0.0,0.0),(0.0,4.0,0.0),(0.0,0.0,3.0)] ::

   begin lattice
   pegloop1 (6.0,4.0) 0.0 0.0
   pegloop1 0.0 (6.0,4.0) 0.0
   indeploop 0.0 0.0 (2.0,3.0)
   end

