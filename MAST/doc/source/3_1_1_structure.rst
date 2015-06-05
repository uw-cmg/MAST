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


