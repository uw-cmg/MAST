Introduction
============
Welcome to the MAterials Simulation Toolkit (MAST)!

**link to whitepaper here??**
MAST is intended to be an easy-to-use wrapper around complex sequences of calculations.

==================
Sequence of events
==================
This is the sequence of events in the MAST program:

#. The user generates an input file, for example, ``test.inp``. See :doc:`Input File <inputfile>`
#. The command ``mast -i test.inp`` uses MAST to parse the input file. MAST then generates a logic tree of ingredients, based on the recipe specified in the inpput file. See :doc:`Recipe <recipe>` and :doc:`Ingredients <ingredients>`.
#. MAST creates a system_recipe_timestamp directory for the particular system and recipe under $MAST_SCRATCH. Each ingredient gets its own directory within the system_recipe_timestamp directory. Recipe information is stored for the scheduling arm of MAST to use.
#. The MAST scheduling arm is run separately, e.g. from a timer. See :doc:`Scheduler <scheduler>`. It runs the ingredients in logical order. When all ingredients in a recipe are complete, the system_recipe_timestamp directory is moved into a $MAST_SCRATCH/complete directory.

==================
Detailed reference
==================

.. toctree::
    :maxdepth: 2
    
    inputfile
    recipe
    ingredients
    scheduler
