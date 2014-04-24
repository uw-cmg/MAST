.. module:: gromacs

=======
Gromacs
=======

Introduction
============

Gromacs is a free classical molecular dynamics package. It is mainly 
used in modeling of biological systems. It is part of the 
ubuntu-linux distribution.
http://www.gromacs.org/

.. warning:: 1) Ase-Gromacs calculator works properly only with gromacs version 4.5.6 or a newer one. (fixed bug related to .g96 format)

.. warning:: 2) It only makes sense to use ase-gromacs for qm/mm or for testing. For pure MM production runs the native gromacs is much much faster (at the moment ase-gromacs has formatted io using .g96 format which is slow).

Gromacs Calculator
==================
This ASE-interface is a preliminary one and it is VERY SLOW so 
do not use it for production runs. It is here because of 
we'll have a QM/MM calculator which is using gromacs as the 
MM part.

For example: (setting for the MM part of a QM/MM run, 
parameter '-nt 1' for serial run)::

  CALC_MM = Gromacs(
    init_structure_file = infile_name,
    structure_file = 'gromacs_qm.g96', \
    force_field='oplsaa', 
    water_model='tip3p',
    base_filename = 'gromacs_qm',
    doing_qmmm = True, freeze_qm = False,
    index_filename = 'index.ndx',
    extra_mdrun_parameters = ' -nt 1 ',
    define = '-DFLEXIBLE',
    integrator = 'md',
    nsteps = '0',
    nstfout = '1',
    nstlog = '1',
    nstenergy = '1',
    nstlist = '1',
    ns_type = 'grid',
    pbc = 'xyz',
    rlist = '1.15',
    coulombtype = 'PME-Switch',
    rcoulomb = '0.8',
    vdwtype = 'shift',
    rvdw = '0.8',
    rvdw_switch = '0.75',
    DispCorr = 'Ener')

For example: (setting for a MM calculation, useful when keeping QM fixed 
and relaxing MM only, parameter '-nt 1' for serial run)::
         
  CALC_MM_RELAX = Gromacs(
    init_structure_file = infile_name,
    structure_file = 'gromacs_mm-relax.g96',
    force_field='oplsaa', 
    water_model='tip3p',
    base_filename = 'gromacs_mm-relax',
    doing_qmmm = False, freeze_qm = True,
    index_filename = 'index.ndx',
    extra_mdrun_parameters = ' -nt 1 ',
    define = '-DFLEXIBLE',
    integrator = 'cg',
    nsteps = '10000',
    nstfout = '10',
    nstlog = '10',
    nstenergy = '10',
    nstlist = '10',
    ns_type = 'grid',
    pbc = 'xyz',
    rlist = '1.15',
    coulombtype = 'PME-Switch',
    rcoulomb = '0.8',
    vdwtype = 'shift',
    rvdw = '0.8',
    rvdw_switch = '0.75',
    DispCorr = 'Ener')


Parameters
==========
The description of the parameters can be found in the Gromacs manual:
http://www.gromacs.org/Documentation/Manual

and extra (ie. non-gromacs) parameter: 

init_structure_file: str
    Name of the input structure file for gromacs
    (only pdb2gmx uses this)
structure_file: str
    Name of the structure file for gromacs
    (in all other context that the initial input file)
force_field: str
    Name of the force field for gromacs
water_model: str
    Name of the water model for gromacs
base_filename: str
    The generated Gromacs file names have this string  
    as the common part in their names (except structure files).
doing_qmmm: logical
    If true we run only single step of gromacs 
    (to get MM forces and energies in QM/MM)
freeze_qm: logical
    If true, the qm atoms will be kept fixed
    (The list of qm atoms is taken from file 'index_filename', below)
index_filename: string
    Name of the index file for gromacs
extra_pdb2gmx_parameters: str
    extra parameter(s) to be passed to gromacs programm 'pdb2gmx'
extra_grompp_parameters: str
    extra parameter(s) to be passed to gromacs program 'grompp'
extra_mdrun_parameters: str
    extra parameter(s) to be passed to gromacs program 'mdrun'

Environmental variables:
========================
  - GMXCMD the name of the main gromacs executable (usually 'mdrun').
    If GMXCMD is not set gromacs test is not run, but in the calculator 
    works using 'mdrun'.
  - GMXCMD_PREF prefix for all gromacs commands (default '')
  - GMXCMD_POST postfix (ie suffix) for all gromacs commands (default '') 
    
Example: MM-only geometry optimization of a histidine molecule
==============================================================
THIS IS NOT A PROPER WAY TO SETUP YOUR MD SIMULATION.
THIS IS JUST A DEMO THAT DOES NOT CRASH. 
(the box size should be iterated by md-NTP, total charge should be 0).

Initial pdb coordinates (file his.pdb):

.. literalinclude:: his.pdb

First generate the initial structure in gromacs format (.gro)

>>> pdb2gmx -f his.pdb -o hish.gro -ff oplsaa -water tip3p 

(pdb2gmx seems strangely to give the molecule hish.gro a total charge of 
-0.110 but this is not a problem for us now, this is just an example)

Then setup a periodic simulation box
 
>>> editconf -f hish.gro -o hish_box.gro -box 3 3 3

Solvate histidine in a water box

>>> genbox -cp hish_box.gro -cs spc216.gro -o hish_waterbox.gro

Generate index file for gromacs groups

>>> make_ndx -f hish_waterbox.gro
>>> q<ENTER>

Finally, relax the structure.
The sample file for relaxation:

.. literalinclude:: gromacs_example_mm_relax.py
