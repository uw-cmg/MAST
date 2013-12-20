===============
The Ingredients
===============
Each ingredient is a separate calculation. Ingredients make up recipes.

Each ingredient is responsible for updating its child ingredients through an ``update_children`` method.

Each ingredient is given:

* A name, which is the full path to the ingredient's directory and is automatically generated using the system name and the recipe template.
* The name of the program running (e.g. vasp)
* A dictionary of program_keys, which contain all program-specific keywords (see :ref:`ingredient-keys`) and come from each ingredient's section in the ``$ingredients`` section of the :doc:`input file <inputfile>`.
* A dictionary of the names of any child ingredients.
* A pymatgen structure object representing the very first structure created from the ``$structure`` section in the input file.

.. _ingredient-keys:

===================
Ingredient Keywords
===================
VASP keywords such as IBRION, ISIF, and so on, can be specified under each ingredient in the ``$ingredients`` section of the :doc:`input file <inputfile>`.

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
* mast_walltime: the walltime requested, in whole number of hours
* mast_memory: the memory per processor requested.
