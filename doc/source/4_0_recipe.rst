#################
The Recipe
#################

********************************
Introduction to the Recipe
********************************

The recipe defines the relationships between ingredients, or which ingredients need to be run before which other ingredients.

Out-of-the-box recipes are stored in ``$MAST_INSTALL_PATH/recipe_templates``. You may copy them into your ``$MAST_RECIPE_PATH`` directory (see :doc:`Installation <1_0_installation>`). If you create new recipes, they should also go in the  ``$MAST_RECIPE_PATH`` directory

The full recipe name goes in the ``$recipe`` section of the input file::

    $recipe
    recipe_file neb.txt
$end

************************************
The Recipe Template
************************************

**Important: when creating or editing recipes, do not use the Tab key. Instead, use 4 spaces to indent.** 

Also make sure that the recipe you are working with has not somehow been converted to tabs.

If you use vi as your code editor, consider adding the following settings to your ``~/.vimrc file``, in order to use python four-space tab stops instead of the Tab character.::

    set tabstop=4
    set shiftwidth=4
    set smarttab
    set expandtab
    set softtabstop=4
    set autoindent


=====================
Syntax
=====================
Each indentation level marks a parent-child relationship.::

    perfect_opt1 (volrelax_lowmesh)
        perfect_opt2
            perfect_opt3
   
The ingredient type of an ingredient is specified in parentheses after the ingredient.

The ingredient type should correspond to ingredient subsections within the $ingredients section of the :doc:`input file <3_0_inputfile>`. If no ingredient type is specified, the ingredient gets all default values from the ingredients_global subsection.

In the recipe::
    
    perfect_opt1 (volrelax_lowmesh)

In the input file::

    $ingredients
    
    begin volrelax_lowmesh
    mast_run_method run_singlerun
    ...
    end
    
    $end


If the parent needs to update several children in different ways, create new trees where the originating parent is the same parent name, but with a different ingredient type. 
*  Those different ingredient types should have different mast_update_children_method keyword values in the input file. 
*  Only the first ingredient type specified per parent, going from the top of the file to the bottom of the file, will be used for all program keywords (run method, write method, INCAR settings, etc.) except for mast_update_children_method. The mast_update_children_method will be taken from the ingredient type specified between the parent and that child. ::

    perfect_stat (stat_to_defect)
        defect_opt
    perfect_stat (stat_to_phonon)
        phonon_opt1

If two children need to be the parent of one ingredient, also create a new tree::

    perfect_stat
        defect_1_opt
        defect_2_opt
    defect_1_opt, defect_2_opt
        neb_1-2_opt

Parent-child relationships are name-based, and the name must also include correct formats for defect labels (defect_XXX), charge labels (q=XX), neb labels (neb_XXX-XXX), and phonon labels (phonon_XXX). These names are important for following the tree structure and for setting the metadata file. Parent-child relationships are specified by these particular folder names. However, once all runs have been completed, post-processing utilities should only look at the metadata file within each run folder, and not at the folder name.

For defects, the labels must correspond to labels in the ``$defects`` section::

    defect_<label>

Defect charges are given as q=p0 for no charge, q=nX for negative charge X (remember that negative charge means more electrons), and q=pX for positive charge X.
(Please note that inducedefect ingredients should be labeled with inducedefect rather than with induce_defect, which will confuse them with defect ingredient labels.)

For nebs, the labels must correspond to labels in the ``$neb`` section::
    
    neb_<label>

For phonons, the labels must correspond to labels in the ``$phonon`` section::

    phonon_<label>
    phonon_<label>_parse
You may create a fully-specified recipe in which you write out the labels, and also the charges, if necessary, for example::
    
defect_opt1_q=n2 (lowmesh)

However, in many cases it is more convenient to use abbreviations within the recipe.
``{begin}`` and ``{end}`` tags specify sections that can be looped over for as many defect labels ``<N>`` are specified in the ``$defects`` section of the input file and NEB labels ``<B-E>``, where ``<B>`` and ``<E>`` are also defect labels, as specified in the ``$neb`` section of the input file.

Charges ``<Q>`` are given by the charge range in the ``$defects`` section. Available charges are carried into the ``<B-E>_<Q>`` labels based on which charges are available to both the ``<B>`` and the ``<E>`` defect in the label. 

Note that defect endpoints need to be the parents of all NEB optimizations and NEB static calculations. 

Example::

    Recipe NEBtest
    perfect_opt1 (lowmesh)
        perfect_opt2
            perfect_stat (static)
            {begin}
            inducedefect_<N> (inducedefect)
                defect_<N>_<Q>_opt1 (lowmesh_defect)
                    defect_<N>_<Q>_opt2 (defect_relax)
                        defect_<N>_<Q>_stat (static)
            {end}
    {begin}
    defect_<N>_<Q>_stat (static)
        phonon_<N>_<Q>_<P> (phonon)
            phonon_<N>_<Q>_<P>_parse (phononparse)
    {end}
    {begin}
    defect_<B>_<Q>_stat (static_to_neb), defect_<E>_<Q>_stat (static_to_neb)
        neb_<B-E>_<Q>_opt1 (neb_to_neb)
            neb_<B-E>_<Q>_opt2 (neb_to_nebstat)
                neb_<B-E>_<Q>_stat (nebstat_to_phonon)
        neb_<B-E>_<Q>_opt2 (neb_to_nebstat)
        neb_<B-E>_<Q>_stat (nebstat_to_phonon)
    {end}
    {begin}
    neb_<B-E>_<Q>_stat (nebstat_to_phonon)
        phonon_<B-E>_<Q>_<P> (phonon)
            phonon_<B-E>_<Q>_<P>_parse (phononparse)
    {end}

