#####################
Programming for MAST
#####################

.. toctree::

   12_0_1_developmentworkflow

===============================
Detailed development workflow
===============================

If you are a MAST programmer in the UW-Madison Computational Materials Group, please see :doc:`12_0_1_developmentworkflow`.

=============================
Source code
=============================

The MAST github repository is located at `https://github.com/uw-cmg/MAST <https://github.com/uw-cmg/MAST>`_.

To report any issues, please create an issue in this repository.

To program with MAST:

#. Clone from the dev branch (see github's instructions for cloning) OR get the latest stable release from `https://github.com/uw-cmg/MAST/releases` and unzip it.

#. Prepend the clone directory (the directory which contains the directory named MAST) to your ``$PYTHONPATH`` environment variable, and the clone directory's ``MAST/bin`` directory to your ``$PATH`` environment variable.

    *  The command ``mast`` should reveal the clone directory instead of any other MAST installation directories.

#. If you previously did a ``pip`` install in order to get the MAST dependencies, go to your python installation's ``lib/python2.7/site-packages/`` directory and rename the MAST package directory there to something else. 
   
    *  Open python and ``import MAST`` at the prompt. ``help(MAST)`` should then reveal MAST in the correct cloned directory, rather than the pip-installed MAST.

To run unit tests and verify that the MAST code is sound, go to the test directory in ``<clone directory>/MAST/test`` and run the command:: 

    nosetests -v --nocapture

The ``nocapture`` option allows print statements.
The ``verbose`` option gives verbose results.

The development team may have designated some tests to be skipped. However, any errors should be reported to the development team as a github issue.

================================
Object hierarchy
================================

Several objects are created in MAST. The classes for these objects are in similarly named files, for example, class MyClass in file myclass.py.

* When the user types ``mast`` or when crontab executes ``mast``, a **MAST monitor** object is created (class MASTmon in MAST). This monitor is responsible for looking through the ``$MAST_SCRATCH`` directory for recipe folders.

* For each recipe folder, 

    * An **Input Options** object is created from the ``input.inp`` file (class InputOptions in MAST/utility, parsed from the input file through class InputParser in MAST/parsers)

    * A **Recipe Plan** object is created from that Input Options object

* The status of the ingredients in the recipe is given by ``status.txt``

    * Depending on the ingredient status, an **Ingredient** object is created using information from the Recipe Plan object (class ChopIngredient, inheriting from class BaseIngredient, in MAST/ingredients)

    * That Ingredient object may involve several **Checker** objects for different programs based on the ``mast_program`` keyword of its ingredient type in :doc:`3_1_2_ingredients` (class XXXChecker, MAST/ingredients/checker)


================================
Code hooks in the input file
================================

The most common modifications to MAST are expected to be:

* Adding support for new programs, e.g. besides VASP

* Adding new parent-child information transfer methods, for example:

    * Giving additional information to a child ingredient, like number of pairs
    
    * Accommodating different run structures, for example, forward on the least symmetric structure among several folders in the parent ingredient

Both of these modifications are currently coded in ``MAST/ingredients/chopingredient.py`` and in ``MAST/ingredients/checker``

In the input file, the ``mast_xxxx_method`` keywords are direct hooks to methods in the **ChopIngredient** class. 

* Methods are separated by semicolons, and can include arguments (see :doc:`3_1_2_ingredients`.)

* The method in the ChopIngredient class may involve a checker, if they are generic but require program-specific treatment, for example, ``forward_final_structure``.

* Or, the method in the ChopIngredient class may not need a checker, if it is totally generic, for example, ``copy_file OLDNAME NEWNAME``

* When used as an update method, please remember that the last argument to a method is going to be the child ingredient's directory, as determined by :doc:`3_1_4_personalrecipe` in the recipe folder.

Support for using a new checker type as self.checker in a ChopIngredient class would need to be added at the top of ``MAST/ingredients/baseingredient.py``.
Alternately, a new checker instance may be initialized on-the-fly within a method, e.g. mychecker = VASPChecker(name=mydirectory)

=========================
Debugging
=========================

For classes which have a self.logger attribute, or functions in which a logger is defined, messages may be logged to the DEBUG level. (self.logger.debug("message"))

Set the MAST_DEBUG environment variable to any value so that the mast.log file will print debug messages.

.. raw:: html

    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-54660326-1', 'auto');
      ga('send', 'pageview');

    </script>

