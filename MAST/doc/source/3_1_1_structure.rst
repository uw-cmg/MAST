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

===================================
Structure indexing (beta)
===================================

Some atomic positions in structures may change significantly after structural relaxation.

Without structure indexing, current operation of MAST is as follows: just before a defect, NEB, or phonon is about to be run, MAST decides which atoms to remove, single out for phonons, etc., based on a coordinate which is guessed at ahead of time, for example, 0.5 0.5 0.5, and sometimes a tolerance.

However, the relaxed atom may be quite far from that coordinate, so the search may fail.

When set to True, the **use_structure_index** keyword in the ``$structure`` section turns on atomic position indexing in the structure.::

    $structure
    use_structure_index True
    ...
    $end

The coordinates in the ``$defect``, ``$neb``, etc. sections should exactly match the coordinates in the initial structure given to MAST. Scaling sizes will be handled automatically.

MAST will create a separate structure index file, a "manifest" file, for the initial structure and each scaling size and defect configuration. Each NEB endpoint will have a separate manifest file, and each phonon calculation will have a separate manifest for detailing selective dynamics information. 

These files are stored in the ``structure_index_files`` directory within the recipe directory.

The manifest file consists of atomic index numbers, which correspond to ``atom_index_<index number>`` files. Each atomic index file is updated with the appropriate coordinates from a completed calculation. 

When a calculation is complete, its relaxed atomic coordinates will be saved back to the appropriate atom_index file.

----------------------------------
NEB notes for structure indexing
----------------------------------
Special care should be taken when defining NEB calculations in the ``$neb`` section of the input file.

Again, vacancy and substitutional defects should have coordinates that exactly match those in the initial structure.

Coordinates in the ``$neb`` section should exactly match coordinates in the ``$defects`` section.

Defect manifests have their defects sorted to the bottom, except for vacancies.
NEB manifests are created from defect manifests.

Atoms indicated in the ``$neb`` section are pulled out of their order in the manifest and put in order at the bottom.

As long as vacancies are accounted for in the ``$neb`` section, and grouped defects are entered in the same order in the defect endpoints in the ``$defects`` section, this process ensures lineup.

Grouped defects that are intended to persist but not move should be placed in the same order, that is, if defect1 and defect2 are going to be in an NEB and both have unmoving Sr antisite defects, then the Sr antisite defects should be in the same order in the begin defect1 and begin defect2 sections.

Pymatgen's interpolate function resorts atoms, so a temporary manifest is made for pymatgen, and the interpolation is resorted back to match the NEB endpoint manifests.

Also, all image atoms are ordered following the endpoint1 manifest. The final endpoint atoms are ordered following the endpoint2 manifest.

Therefore, after completion, coordinates for image atoms are recorded in the atom index files indicated by the endpoint 1 manifest.


.. raw:: html

    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-54660326-1', 'auto');
      ga('send', 'pageview');

    </script>

