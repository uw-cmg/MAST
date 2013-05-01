Tutorial
========
Welcome to the MAterials Simulation Toolkit (MAST)!

============
Installation
============
#. Copy the top MAST directory into your user directory. We will refer to the top MAST directory as ``//home/user/topmast`` in this documentation.

#. Set the following environment variables. You may want to put the export commands in your setup profile, such as ``//home/user/.bashrc`` and then either run ``source ~/.bashrc`` or log out and log back in to your user terminal.

    #. MAST_INSTALL_PATH: This variable should be set to the installation directory, for example::
    
        export MAST_INSTALL_PATH=//home/user/topmast

    #. MAST_RECIPE_PATH: Initially, this variable should be set to the existing recipe_templates folder. As you develop additional recipes, you may want to change this variable to reflect your customized organization of recipes.::
        
        export MAST_RECIPE_PATH=//home/user/topmast/recipe_templates

    #. MAST_SCRATCH: This variable may be set to any directory. Default behavior is to generate ingredients under this MAST_SCRATCH directory, unless the match_scratch keyword is specified with an overriding path in the input file.

    #. PYTHONPATH: If this environment variable already exists, the installation directory should be appended. Otherwise, this variable can be set to the installation directory. Assuming PYTHONPATH already has some value  (use ``env`` to see a list of environment variables):: 
        
        export PYTHONPATH=$PYTHONPATH://home/user/topmast
        
    #. VASP_PSP_DIR: This variable is necessary if VASP and VASP pseudopotential files are being used. See the documentation for the Materials Project's **pymatgen** code.
    #. PATH: This variable should be appended with the MAST bin directory, for example::
    
        export PATH=$PATH://home/user/topmast/bin

The script ``initialize.py`` in ``//home/user/topmast`` will give a set of default values to paste into your terminal window or ~/.bashrc file.

==================
Sequence of events
==================
This is the sequence of events in the MAST program:

#. The user generates an input file, for example, ``test.inp``. See :ref:`input-file`
#. The command ``mast -i test.inp`` uses MAST to parse the input file. MAST then generates a logic tree of ingredients, based on the recipe specified in the inpput file. See :ref:`recipe` and :ref:`ingredients`.
#. MAST creates a system_recipe_timestamp directory for the particular system and recipe under $MAST_SCRATCH. Each ingredient gets its own directory within the system_recipe_timestamp directory. Recipe information is stored for the scheduling arm of MAST to use.
#. The MAST scheduling arm is run separately, e.g. from a timer. See :ref:`scheduler`. It runs the ingredients in logical order. When all ingredients in a recipe are complete, the system_recipe_timestamp directory is moved into a $MAST_SCRATCH/complete directory.


.. _input-file:

===============
Input File
===============
The input file contains several sections denoted by ``$sectionname`` and ``$end``
Here is a sample input file, ``test.inp``::
    
    # Test file
    $mast
    program VASP
    system_name MgAl
    $end

    $structure
    coord_type fractional
    begin coordinates
    Mg 0.000000 0.000000 0.000000
    Al 0.500000 0.500000 0.500000
    end

    begin lattice
    3.0 0.0 0.0
    0.0 3.0 0.0
    0.0 0.0 3.0
    end
    $end

    $ingredients
    begin ingredients_global
    mast_kpoints 3x3x3
    mast_xc pbe
    end

    begin optimize
    encut 300
    ibrion 2
    end
    $end

    $recipe
    recipe_file recipe_test.txt
    $end

The ``$mast`` section contains the program and system name. An optional mast_scratch keyword and path to a directory may be given to write the recipe and ingredients under this directory rather than under the directory given in $MAST_SCRATCH.

The ``$structure`` section contains the coordinate type, coordinates, and lattice. Optionally, the name of a VASP POSCAR-type file can be inserted here using the keyword posfile, e.g. ``posfile fcc_POSCAR``. The coordinates are given with element name and then three fractional coordinates along the lattice vectors.

The ``$ingredients`` section contains a section for global ingredient keywords and then a section for each separate ingredient. VASP INCAR keywords are included in these sections. All other keywords are prefaced with ``mast_``. A listing of available keywords is in :ref:`ingredient-keys`.

The ``$recipe`` section contains a section for the recipe template to be used.

Other sections include:

* The ``$defects`` section, which includes the defect type (vacancy or interstitial), the defect coordinates, and the defect element symbol::
    
    $defects
    vacancy 0 0 0 Mg
    vacancy 0.5 0.5 0.5 Mg
    interstitial 0.25 0.25 0 Mg
    interstitial 0.25 0.75 0 Mg
    $end

* The ``$neb`` section, which includes a list of nudged-elastic-band hops, corresponding to the defects listed in the ``$defects`` section, and the number of interpolated images for each hop. For example,::

    $neb
    hops 1-2 1-3 3-4
    images 3
    $end

.. _recipe:

=============
The Recipe
=============
Here is an example recipe template::

    Recipe TEST

    Ingredient <sys>_perfect_opt1 Optimize
    Ingredient <sys>_perfect_opt2 Optimize

    Parent <sys>_perfect_opt1 child <sys>_perfect_opt2::structure

The recipe contains:
#. The recipe name
#. Each ingredient in order, including the desired ingredient name format, and the ingredient type
#. Logical relationships between Parent and Children ingredients, and the information passed to the child. <<LAST PART DEPRECATED?>>

* <sys> will be replaced with the system name from the input file.
* <N> will be replaced with defect numbers, in order, with NEB hops, and with NEB image numbers.

.. _ingredients:

===============
The Ingredients
===============
Each ingredient is a separate calculation. Ingredients make up recipes.

Each ingredient is responsible for updating its child ingredients through an ``update_children`` method.

Each ingredient is given:

* A name, which is the full path to the ingredient's directory and is automatically generated using the system name and the recipe template.
* The name of the program running (e.g. vasp)
* A dictionary of program_keys, which contain all program-specific keywords (see :ref:`ingredient-keys`) and come from each ingredient's section in the ``$ingredients`` section of the input file
* A dictionary of the names of any child ingredients.
* A pymatgen structure object representing the very first structure created from the ``$structure`` section in the input file.

.. _scheduler:

===============
The Scheduler
===============
The scheduler ensures that ingredients run in the correct order.

.. _ingredient-keys:

===================
Ingredient Keywords
===================
VASP keywords such as IBRION, ISIF, and so on, can be specified under each ingredient in the ``$ingredients`` section of the input file.

Any keyword prepended by ``mast_`` is considered a special keyword and will not be written into the VASP INCAR.


* mast_kpoints: specify k-point instructions in the form of kpoints along lattice vectors a, b, and c, and then a designation M for Monkhorst-Pack or G for Gamma-centered. ``mast_kpoints = 3x3x3 G``
    * required for VASP

* mast_xc: specify an exchange correlation functional; for VASP, follow the convetions of pymatgen.
    * required for VASP

* mast_multiplyencut: specify a number with which to multiply the maximum ENCUT value of the pseudopotentials. Volume relaxations in VASP take 1.5; otherwise 1.25 is sufficient.
    * defaults to 1.5

* mast_setmagmom: specify a string to use for setting the initial magnetic moment. A short string will result in multipliers, ex: 1 5 1 = 2*1 2*5 8*1 for a 12-atom unit cell. A string of the number of atoms in the POSCAR will be printed as entered, for example, 1 -1 1 -1 1 -1 1 -1.

* mast_adjustnelect: specify an adjustment to the total number of electrons. For example, -2 to remove two electrons, and +2 to add two electrons

The following queue-submission keywords are discussed more in :ref:`platforms`. 

* mast_processors: the total number of processors requested. Use this or use mast_nodes and mast_ppn.
* mast_nodes: the number of nodes requested.
* mast_ppn: the number of processors per node requested.
* mast_queue: the queue requested.
* mast_exec: the full executable line, including any mpi commands. Remember that all input file options are turned into lowercase, with the exception of mast_xc, which is turned into all uppercase in the INCAR file.
* mast_walltime: the walltime requested.
* mast_memory: the memory per processor requested.

.. _platforms:

================
Platform Support
================
Queue and submission script commands are in ``//home/user/topmast/submit`` and may need to be heavily modified depending on the platform used. 
To customize the queue submission behavior, copy the two _example.py files into new .py files, removing "_example" and overwriting the default files.
 
The out-of-the-box PBS submission script is built using

* mast_processors or a combination of mast_ppn and mast_nodes
* mast_queue
* mast_exec
* mast_walltime
* mast_memory
* the ingredient name
