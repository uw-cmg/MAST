#############################
What's new in version 1.2.X
#############################

**Additions:**

* Finite size scaling support (the <S> tag) has been added. See :doc:`3_1_1_structure`, :doc:`3_1_3_recipe`, and :doc:`6_0_postprocessingtools`.

**Fixes:**

* The :doc:`8_0_2_gbdiff` and :doc:`8_0_3_diffanalyzer` packages are now properly included in the installation directory after running ``setup.py``. 

**Changes for users:**

* The ``$recipe`` section of the input file now requires the recipe to be entered directly.

    * Do not use a text file name any more.

    * Do not start with a recipe name line.

    * The ``MAST_RECIPE_PATH`` environment variable is no longer necessary.

* When the input file is processed, it will create a ``$personal_recipe`` section directly in the input file.
    
    * There is no longer a ``personal_recipe.txt`` file in the recipe directory.

    * If copying an input file for use in a new recipe, delete the ``$personal_recipe`` section from the new copy of the input file.

* MAST will now tell you where it was installed when you run ``mast``.

* Platform support is now all under ``<MAST installation directory>/submit/platforms/``. 

    * The ``platforms`` folder is no longer copied to ``$MAST_CONTROL``. See :doc:`1_0_installation` for creating and modifying platforms.


**Changes for programmers:**

* Automatic citation support files are no longer copied to ``$MAST_CONTROL``. They are located in ``<MAST installation directory>/summary/citations``.

* Program key files are no longer copied to ``$MAST_CONTROL``. They are located in ``<MAST installation directory>/ingredients/programkeys``.

* Optimizer.py is no longer copied to ``$MAST_CONTROL``. It is located in ``<MAST installation directory>/structopt``

* mastmon_submit.sh is now copied from ``<MAST installation directory>/submit/platforms/<platform name>`` into ``$MAST_CONTROL`` each time that MAST is run. Edits should therefore be made to ``<MAST installation directory>/submit/platforms/<platform name>/mastmon_submit.sh`` if they are necessary.

* The ``$MAST_CONTROL/set_platform`` file is no longer used and references to it have been removed.
