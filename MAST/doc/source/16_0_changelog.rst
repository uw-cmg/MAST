#############################
Changelog
#############################

==============
version 2.0.1
==============
* Fixes to defect formation energy tool for compatibility with pymatgen 4.7.1
* Pin pymatgen version in setup.py to 4.7.1
* Replace "what's new" files with this changelog

==============
version 2.0.0
==============

**Additions:**

None.

**Fixes:**

* MAST is now compatible with pymatgen 4.X.X
* Major fixes have been made to atom indexing.

**Changes for users:**

* StructOpt packaged with MAST is now deprecated. Please see the code base maintained at `StructOpt_modular <https://github.com/uw-cmg/StructOpt_modular>`_ and/or contact its developers for assistance.

Users with old workflows combining StructOpt and MAST should use version 1.3.4.

**Changes for programmers:**

None.

=========================
version 1.3.0
=========================

**Additions:**

* STEM image simulation functionality has been added to structopt. See :doc:`8_0_4_structopt`.
* A charge-density-based pathfinding method has been added as an alternative to linear interpolation for NEBs. See :doc:`3_1_2_ingredients`.
* An optional atom indexing feature has been implemented for structures which are expected to undergo significant atomic relaxation when relaxed and/or defected. See :doc:`3_1_1_structure`.
* A 14-frequency concentrated diffusion model has been added to the Diffusion Coefficient post-processing tool. See :doc:`6_0_postprocessingtools`.
* The new MAST "modify recipe" function allows existing recipes to be modified. See :doc:`5_0_runningmast`.

**Fixes:**


**Changes for users:**

* Scaling is its own section in the input file, and is no longer under the structure section. See :doc:`3_1_9_scaling` for format changes.
* Instead of a mast_recipe.log file in each recipe directory, there is now one large mast.log file in $MAST_CONTROL.
* Instead of recipe-level trapped errors causing MAST to fail, MAST will create a MAST_ERROR document in the recipe directory and log a warning to ``$MAST_CONTROL/mast.log``

* Standalone grain-boundary diffusion and diffusion analyzer tools must be unzipped from the MAST tar.gz file downloaded from github (see :doc:`12_0_programming`). They will not be installed through pip.

**Changes for programmers:**

* In order to see some logged messages that are now logged to the DEBUG level, set environment variable MAST_DEBUG to any value.

====================
version 1.2
====================

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

.. raw:: html

    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-54660326-1', 'auto');
      ga('send', 'pageview');

    </script>

