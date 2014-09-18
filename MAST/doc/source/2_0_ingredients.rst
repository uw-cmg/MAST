########################
Ingredients
########################

Each ingredient is a separate calculation. Ingredients make up recipes.

Each ingredient is responsible for updating its child ingredients through an update_children method.

Each ingredient is given:

*  A name, which is the full path to the ingredient's directory and is automatically generated from information in the input file.

    *  Some ingredient names must be structured specifically. For examples of naming conventions, see the :doc:`Recipe  <4_0_recipe>`. In particular:

    *  An ingredient which is supposed to correspond to values given by the ``$defects`` section of the :doc:`Input File <3_0_inputfile>` should always be named with ``inducedefect_`` (for the structural creation of the defect) or ``defect_`` (for an actual defect calculation)

    *  An ingredient which is supposed to correspond to values in the ``$neb`` section, such as a nudged elastic band (NEB) calculation or the static image calculations of an NEB calculation, should always be named with ``neb_``

    *  A phonon calculation should always be named with ``phonon_``, and a subsequent calculation of phonon results should be named with ``phonon_...parse``

    *  The letters **q=** are reserved (generated automatically by the recipe template in some cases) and should not otherwise be put in an ingredient name
    

*  A dictionary of program-specific keywords, which come from each ingredient.s section in the ``$ingredients`` section of the :doc:`Input File <3_0_inputfile>`.

*  A pymatgen structure object representing the very first structure created from the ``$structure`` section in the input file.

*  A type, which is specified in the recipe, next to the ingredient name, in parentheses. The ingredient type corresponds to the ingredient type subsection in the ``$ingredients`` section of the input file. The information given in these subsections includes:

    *  Program-specific keywords

    *  Other MAST keywords, including:

        *  The **write** method: which files the ingredient should write out before running (e.g., create the INCAR)

        *  The **ready** method: how MAST can tell if the ingredient is ready to run (often, in addition to writing its own files, an ingredient must also wait for data from its parent ingredient(s)). 
            
        *  The **run** method: what MAST should do to run the ingredient (e.g. submit a submission script to a queue, or perform some other action)
            
        *  The **complete** method: how MAST can tell if the ingredient is considered complete
            
        *  The **update children** method: what information an ingredient passes on to its children, and how this information is passed on

The same ingredient in a recipe may be listed more than once, with several different ingredient types. In this case, the first four methods and all the ingredient keywords are given by the first ingredient type encountered. Only the .update_children. method is changed for all subsequent positions. This situation indicates that the ingredient has many children, which must be updated in different ways and thus needs different update_children methods for those different situations.

More detail on ingredients is given in the ``$ingredients`` section of the :doc:`Input File <3_0_inputfile>`.

