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
    * A **Recipe Plan** object is created from that Input Options object and from the ``personal_recipe.txt`` file (class RecipePlan in MAST/recipe)

* The status of the ingredients in the recipe is given by ``status.txt``

    * Depending on the ingredient status, an **Ingredient** object is created using information from the Recipe Plan object (class ChopIngredient, inheriting from class BaseIngredient, in MAST/ingredients)
    * That Ingredient object may involve several **Checker** objects for different programs (class XXXChecker, MAST/ingredients/checker)


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

* Methods are separated by semicolons, and can include arguments (see :ref:`3_0_inputfile`)
* The method in the ChopIngredient class may involve a checker, if they are generic but require program-specific treatment, for example, ``forward_final_structure``.
* Or, the method in the ChopIngredient class may not need a checker, if it is totally generic, for example, ``copy_file OLDNAME NEWNAME``
* When used as an update method, please remember that the last argument to a method is going to be the child ingredient's directory, as determined by ``personal_recipe.txt`` in the recipe folder.

Support for using a new checker type as self.checker in a ChopIngredient class would need to be added at the top of ``MAST/ingredients/baseingredient.py``.
Alternately, a new checker instance may be initialized on-the-fly within a method, e.g. mychecker = VASPChecker(name=mydirectory)

--------------------------------
Adding a new checker
--------------------------------
This section contains information pertinant to creating a new checker connecting a program with MAST

* What is a Checker

    * Object class that contains all of the functions relative to setting up, running, submitting, and evaluating output for a program

    * Essentially, the checker is a python wrapper for communicating with a specific program

* What are the main functions and keywords in a checker

    * all checkers should have the same basic initialization function

        * should be initialized with the name, program_keys, and structure keywords
            * name is a string which contains the name of the directory in which the program will be run
                * this will correspond to a specific ingredient from the MAST input file
            * program_keys is a dictionary which contains all of the keyword from the ingredient section of the MAST input file
            * structure is a pymatgen atomic structure object that is inherited from the MAST input file to be used by the program
        * the checker is also initialized with the BaseChecker class
            * this class contains many standard functions relevant to a program
                * notably including: writing submit scripts, running the program, changing the status, etc.
            * note that you will need to update the __init__ function in the BaseChecker class with your new checker and program or else rely on the generic checker

    * all checkers should have the following functions in order to integrate with MAST:
        * is_complete -> checks to see if program has finished running
        * is_frozen -> checks to see if program is froze
        * is_ready_to_run - -> checks to see if necessary files are ready for program to run
        * is_started -> checks to see if program has started running
        * set_up_program_input -> sets up the required input files for the program
        * forward_final_structure_from_directory -> copies the output file from one directory to the input of another
        * get_structure_from_file -> gets a structure from a file related to the program
        * _get_allowed_keywords -> checks what keywords are specific for the program
        * _get_non_mast_keywords -> gets only the program specific keywords from the MAST input file

    * adding a new checker also requires adding a program_allowed_keywords.py file to the MAST/ingredients/programkeys folder
        * should include only the keywords allowed in the program
        * need to be careful that these don't overlap other programs that might be used by other programs called by the same MAST input file

* How the checker gets called (How to make the program go)

    * the checker is called by the chopingredient.py script
        * chop ingredient handles all of the specifics for a calculation
        * the methods used by the chop ingredient script are set by the mast commands in the ingredients section of the MAST Input file
        