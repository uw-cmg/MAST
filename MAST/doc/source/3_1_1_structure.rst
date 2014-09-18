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

The file should be located in the same directory as the input file at the time you call MAST, and should not be moved until the recipe is complete.

A CIF file should end with .cif.

A POSCAR-type filename must start with ``POSCAR_`` or ``CONTCAR_`` in order for pymatgen to recognize it. The elements will be obtained from the POSCAR unless you also have a POTCAR in the directory, in which case, check your output carefully because the elements might be given by the POTCAR instead, no matter what elements are written in the POSCAR file.

====================================
Structure by specification
====================================

To specify a structure, use the following subsections:

**coord_type**: This keyword specifies fractional or cartesian coordinates. Only fractional coordinates have been thoroughly tested with most MAST features.

**lattice**: The lattice subsection specifies lattice basis vectors on a cartesian coordinate system.

**elementmap**: The elementmap subsection allows you to create a generic lattice and interchange other elements onto it. This is useful when looping over other elements (discussed in :doc:`3_0_inputfile`).

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


================================
Finite-size scaling
================================
Finite size scaling is supported with a special "scaling" subsection.

Defect positions will be automatically scaled.

    * For example, ``0.25 0.0 0.0`` in the original supercell would become ``0.125 0.0 0.0`` in a 2x1x1 cell. 

Special notes:

*  :doc:`3_1_2_ingredients` should include an "inducescaling" ingredient with a ``mast_run_method`` of ``run_scale``

*  :doc:`3_1_3_recipe` should include ``inducescaling_<S>`` and ``defect_<S>`` ingredients.

    *  The "<S>" tags will correspond to the scaling sizes and labels.

The scaling section has the syntax:

    * Scaling matrix of integers ``[M, N, P]`` or ``[M1 M2 M3, N1 N2 N3, P1 P2 P3]``
    
    * Kpoint mesh in the form ``QxRxS``

    * Kpoint mesh type, M for Monkhorst-Pack and G for Gamma-point centered

    * (Optional) Kpoint mesh shift, in floats, e.g. ``0.1 0.2 0.3``

    * (Optional) Label, in the form ``label=<labelname>``

Example::
  
    begin scaling
    [1 0 0,0 1 0, 0 0 1] 4x4x4 M label=1x1x1
    [2 0 0,0 2 0, 0 0 1] 2x2x4 M label=2x2x1
    [2 0 0,0 2 0, 0 0 2] 2x2x2 M label=2x2x2
    [3 0 0,0 3 0, 0 0 3] 1x1x1 M label=3x3x3
    end

In order to figure out which scaling sizes to use for finite-size scaling, MAST includes a Madelung potential utility.

This utility generates a distribution of cell sizes for best scaling, according to the method in::

    Hine, N. D. M., Frensch, K., Foulkes, W. M. C. & Finnis, M. W. Supercell size scaling of density functional theory formation energies of charged defects. Physical Review B 79, 13, doi:10.1103/PhysRevB.79.024112 (2009).

Run this utility as follows in order to generate a cut-and-paste for the scaling section. ::

    mast_finite_size_scaling_sizes perfDir defDir minDefDist maxNumAtoms numStructAsked

* **perfDir**: perfect primordial (small) cell directory, which should already have run and include VASP CONTCAR, OSZICAR, etc. files.

* **defDir**: defected primordial cell directory, which should already have run and include VASP CONTCAR, OSZICAR, etc. files.

* **minDefDist** (default 3): minimum defect-defect distance between periodic images, in Angstroms.

* **maxNumAtoms** (default 600): maximum number of atoms for supercell size evaluations

* **numStructAsked** (default 5): number of structures to return in the distribution 

* Note that you will have to manually adjust the kpoint mesh in your cut-and-paste.

