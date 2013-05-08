=============
The Recipe
=============
Here is an example recipe template::

    Recipe TEST

    Ingredient <sys>_perfect_opt1 Optimize
    Ingredient <sys>_perfect_opt2 Optimize

    Parent <sys>_perfect_opt1 child <sys>_perfect_opt2::structure

The recipe contains:
#. The recipe name
#. Each ingredient in order, including the desired ingredient name format, and the ingredient type
#. Logical relationships between Parent and Children ingredients, and the information passed to the child. <<LAST PART DEPRECATED?>>  A parent ingredient may have more than one child, and a child ingredient may have more than one parent. Group multiple parents and children under one keyword, for example, ``Parent A B Child C``.

* <sys> will be replaced with the system name from the input file.
* <N> will be replaced with defect numbers, in order, with NEB hops, and with NEB image numbers.
