Introduction
============
Welcome to the MAterials Simulation Toolkit (MAST)!

MAST is an automated workflow manager and post-processing tool. 

MAST focuses on diffusion and defect workflows that use density functional theory. It interfaces primarily with the Vienna Ab-initio Simulation Package (VASP). 

However, MAST can be generalized to other workflows and codes.

`MAST is available from the Python Package Index. <https://pypi.python.org/pypi/MAST>`_

Additional tools and unit tests are available through the latest `MAST tar.gz file. <https://github.com/uw-cmg/MAST/releases>`_

==================
The MAST Kitchen
==================

MAST uses kitchen terminology to organize the materials simulation workflow.

* An :doc:`Ingredient<2_0_ingredients>` is a single calculation, like a single VASP calculation resulting in a relaxed structure and energy. 

* A :doc:`Recipe<4_0_recipe>` is a collection of several ingredients, including information about how the ingredients are combined together. 

    * As in a cooking recipe, ingredients may need to be addressed in a logical order, with some ingredients depending on other ingredients.

    * :doc:`3_1_3_recipe` defines this order, or workflow.

When MAST reads an input file, it creates a recipe in the ``$MAST_SCRATCH`` directory.

* Many recipes can reside in ``$MAST_SCRATCH``.
    
* MAST will check and update the recipes in alphanumeric order.

When MAST finds that a recipe is complete, it will move the recipe from ``$MAST_SCRATCH`` to ``$MAST_ARCHIVE``.

=============================
Computing in the MAST Kitchen
=============================

#.  Install MAST (see :doc:`1_0_installation`).

#.  Plan your workflow. 

    * What are the single calculations you will need (Ingredients)? 

    * Which calculations depend on each other and should be grouped into a Recipe? 
    * What are all of the conditions for each calculation?
    
        * Which calculations have a volume change?
        
        * Which calculations should be run at fixed volume?
        
        * How fine a kpoint mesh does each calculation need?
        
        * Etc...

#.  Run an example file (see :doc:`17_0_testmast`) to get a feel for how MAST works.

#.  Copy and modify an example file for your own workflow.

**Please check your output carefully, especially when setting up a new workflow.**
