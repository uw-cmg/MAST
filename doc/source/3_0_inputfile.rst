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

The ``$ingredients`` section contains a section for global ingredient keywords and then a section for each ingredient type. VASP INCAR keywords are included in these sections. All other keywords are prefaced with ``mast_``. A listing of available keywords is in the :doc:`Ingredients <ingredients>`. Each ingredient type in the recipe should have a begin ingredienttype, end section, even if there are no keywords within that section.

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
