########################
Recipes
########################

Each recipe is a collection of ingredients.

The recipe directory will contain:

* An ``input.inp`` file, which is a copy of the original input file and is used by MAST when checking the recipe. The original input file is not used. This copy also contains :doc:`3_1_4_personalrecipe`, which is not in the original input file.

* Archive files from the initial setup of the recipe directory

* Ingredient directories

* A top-level ``metadata.txt`` file, which stores important information for MAST

* A ``status.txt`` file listing the status of each ingredient

* If an error is detected, a ``MAST_ERROR`` file will be created using the error text.

For other logging information, see the ``$MAST_CONTROL/mast.log`` file. 

A recipe object (created by MAST from the input file, and accessible to MAST while MAST is running) will have:

*  A name, which is the full path to the recipe's directory

*  Several dictionaries which specify:

    *  Which ingredient directories exist
    
    *  Which ingredients have parents, and the names of those parent ingredients
    
    *  Which method(s) each ingredient should run for each mast_xxx_method (see :doc:`3_1_2_ingredients`)

        * Which method(s) each ingredient should run for its mast_update_children_method, depending on the name of the child ingredient

