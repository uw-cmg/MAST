Introduction
============
Welcome to the MAterials Simulation Toolkit (MAST)!

**link to whitepaper here??**
MAST is intended to be an easy-to-use wrapper around complex sequences of calculations.

==================
The MAST Kitchen
==================

MAST uses kitchen terminology to organize the materials simulation workflow.

* An :doc:`Ingredient <ingredients>` is a single calculation, like a single VASP calculation resulting in a relaxed structure and energy. 

* A :doc:`Recipe <recipe>` is a collection of several ingredients. As in a cooking recipe, ingredients may need to be addressed in a certain logical order. The recipe template :doc:`Recipe Template <recipetemplate>` defines this order.

* A :doc:`Buffet <buffet>` is a collection of recipes. The buffet groups recipes which depend on one another. Like in a buffet line, MAST will address recipes one after another.

=============================
Computing in the MAST Kitchen
=============================

#. Plan your workflow. What are the single calculations you will need (Ingredients)? Which calculations depend on each other and should be grouped into a Recipe? Will all of your Recipes be independent, or do some Recipes depend on the output of other Recipes?

#. Create a new recipe or use some of the standard recipes in your $MAST_RECIPE_TEMPLATES directory. See :doc:`Recipe Template <recipetemplate>` for help.

#. Create an input file, for example, ``test.inp`` See :doc:`Input File <inputfile>`

#. Run the command ``mast -i test.inp`` to parse the input file. Under $MAST_SCRATCH, MAST creates a buffet_timestamp directory, and within that directory, a system_recipe_timestamp directory for the particular system and recipe. Each ingredient gets its own directory within the system_recipe_timestamp directory. The buffet information, including all of the input options, is stored in a :doc:`Buffet Pickle <buffetpickle>`.

#. Run the command ``mast`` to start the MAST scheduling arm. See :doc:`Scheduler <scheduler>`. The MAST scheduler will get information from the buffet pickles, look at each recipe in the buffet, and run the recipe's ingredients. When all ingredients in all recipes of the buffet are complete, the buffet_timestamp directory is moved into a $MAST_SCRATCH/archive directory.

**link to a diagram of the MAST workflow**

==================
Detailed reference
==================

.. toctree::
    :maxdepth: 2
    
    inputfile
    recipe
    ingredients
    scheduler
    repetitions
