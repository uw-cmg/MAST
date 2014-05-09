Introduction
============
Welcome to the MAterials Simulation Toolkit (MAST)!

MAST is intended to be an easy-to-use wrapper to facilitate complex sequences of calculations.

==================
The MAST Kitchen
==================

MAST uses kitchen terminology to organize the materials simulation workflow.

* An :doc:`Ingredient <2_0_ingredients>` is a single calculation, like a single VASP calculation resulting in a relaxed structure and energy. 
* A :doc:`Recipe <4_0_recipe>` is a collection of several ingredients and information about how the ingredients are combined together. 

    * As in a cooking recipe, ingredients may need to be addressed in a certain logical order. This temporal order of how ingredients work together is the workflow.
    * The :doc:`Recipe Template <4_0_recipe>` and :doc:`Input File <3_0_inputfile>` together describe the order of the ingredients and the way they are combined together. 

=============================
Computing in the MAST Kitchen
=============================

#.  Install MAST (see :ref:`1_0_installation`).

#.  Plan your workflow. 

    * What are the single calculations you will need (Ingredients)? 
    * Which calculations depend on each other and should be grouped into a Recipe? 
    * What are all of the conditions for each calculation (e.g. which ones can have volume change, and which ones should be at fixed volume? How fine a kpoint mesh does each calculation need? etc.)?

#.  Start with some of the standard recipes in your $MAST_RECIPE_PATH directory or use a new template.
#.  Create an input file, for example, test.inp.
#.  Run the command ``mast -i test.inp`` to parse the input file. 
#.  Under ``$MAST_SCRATCH``, MAST creates a timestamped recipe directory. 
#.  Within the recipe directory:

    #.  Each ingredient gets its own directory within the system_recipe_timestamp directory.
    #.  Additional files are created, including:

        #. ``personal_recipe.txt``, which is your recipe template file filled in with information gathered from the input ``.inp`` file.
        #. ``archive_input_options.txt``, so you can see what the input options originally were
        #. ``archive_recipe_plan.txt``, which tells you how MAST interpreted the recipe file. You can check this file to see which ingredients are considered parents of which other ingredients, for troubleshootin
        #. ``status.txt``, which tells the status of all the ingredients.
        #. ``input.inp``, which is a copy of the input file (or an individual loop of a looped input file)
        #. ``metadata.txt``, which stores metadata information
        #. ``mast_recipe.log``, which stores recipe-level logging information.

#.  Run the command ``mast`` to start the MAST scheduling arm. The MAST scheduler will get information from the personal_recipe.txt, input.inp, and status.txt file in the recipe folder.
#.  When all ingredients in the recipe are complete, the recipe directory is moved into a ``$MAST_ARCHIVE`` directory.

Please check your output carefully, especially when setting up a new workflow using MAST.
