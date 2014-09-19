#####################
Programming for MAST
#####################
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
