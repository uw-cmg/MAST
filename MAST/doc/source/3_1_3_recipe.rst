####################
The Recipe section
####################

The ``$recipe`` section of the input file contains information about how the ingredients are related to each other.

*  This information complements the ``mast_update_children_method`` keyword given for each ingredient.

An ingredient in the recipe is referred to by::

    <ingredient name> (ingredient type in $ingredients section)

For example::

    perfect_opt1 (lowmesh_relaxation)

If no ingredient type is given, then only settings from the ingredients_global ingredient type of the input file will be used.

The ingredient name has some restrictions:

    *  For a simple workflow, the ingredient name may be fully and arbitrarily specified for the user.

    *  In most complex workflows, however, tags may be used as shortcuts to ingredient names. These tags will be filled in from information in the input file.

        * **<S>**: The ``scaling`` subsection of :doc:`3_1_1_structure`
        * **<N>**: :doc:`3_1_5_defects`
        * **<Q>**: The ``charge`` keyword in :doc:`3_1_5_defects`
        * **<P>**: The ``phonon`` keyword in :doc:`3_1_5_defects` and `3_1_6_neb`
        * **<B>, <E>, <B-E>**: :doc:`3_1_6_neb`

        * The filled-in tags will be evident in :doc:`3_1_4_personalrecipe` of the ``input.inp`` file in the recipe directory, once MAST has read the input file and set up the recipe directory.

    *  When tags are used, certain conventions must be followed:
    
        * Inducing scaling must use an ``inducescaling_<S>`` ingredient.

        * Inducing defects must use an ``inducedefect_<N>`` or ``inducedefect_<S>_<N>`` ingredient.

        * Defects must start with ``defect``, and if tags are used, they must follow the order <S>, <N, B, or E>, <Q>, depending on which tags are used. ::
        
            defect_<S>_<N>_<Q>_arbitrarysuffix
        
        * Phonons must start with ``phonon``, and if tags are used, they must follow the order <S>, <N or B-E>, <Q>, <P>

        * NEBs must start with ``neb``, and if tags are used, they must follow the order <S>, <B-E>, <Q>
    
**Important: when creating or editing recipes, do not use the Tab key. Instead, use 4 spaces to indent.** 

    * See :doc:`1_0_installation` for setting up text editors.

    * Also make sure that the recipe you are working with has not somehow been converted to tabs.

=====================
Syntax
=====================
Each indentation level marks a parent-child relationship.::

    perfect_opt1 (volrelax_lowmesh)
        perfect_opt2
            perfect_opt3
   
The ingredient type of an ingredient is specified in parentheses after the ingredient.

The ingredient type should correspond to ingredient subsections within :doc:`3_1_2_ingredients`. If no ingredient type is specified, the ingredient gets all default values from the ingredients_global subsection.

In the recipe::
    
    perfect_opt1 (volrelax_lowmesh)

In the input file::

    $ingredients
    
    begin volrelax_lowmesh
    mast_run_method run_singlerun
    ...
    end
    
    $end


If the parent needs to update several children in different ways, create new trees where the originating parent is the same parent name, but with a different ingredient type::
    
    perfect_stat (stat_to_defect)
        defect_opt
    perfect_stat (stat_to_phonon)
        phonon_opt1


*  Those different ingredient types should have different mast_update_children_method keyword values in the input file. 

*  They should have all the same other keywords.

If two children need to be the parent of one ingredient, also create a new tree::

    perfect_stat
        defect_1_opt
        defect_2_opt
    defect_1_opt, defect_2_opt
        neb_1-2_opt

Parent-child relationships are name-based, and the name must also include correct formats for size-scaling labels <S>, defect labels <N, B, or E>, neb labels <B-E>, charge labels <Q>, and phonon labels <P>.

*  These names are important for following the tree structure and for setting the metadata file. 
*  Parent-child relationships are specified by these particular folder names.
*  Some post-processing utilities may also rely on folder names.

**The <S> tag**
The <S> tag will correspond to labels in the ``scaling`` subsection of :doc:`3_1_1_structure`.

**The <N>, <B>, <E>, and <B-E> tags**
For defects, the <N> tag will correspond to labels in :doc:`3_1_5_defects`.

The same labels will be matched up and should be used as <B> and <E> labels (beginning and ending states) to correspond with NEBs, which are labeled <B-E>.

The NEB labels will correspond to labels in :doc:`3_1_6_neb`

NEB label names must match up exactly with defect label names. For example, defect_vac1 and defect_vac2 must match up with neb_vac1-vac2.

Use <N> in a recipe unless specifying that a defect is a parent of an NEB, in which case use <B> or <E>::

    {begin}
    defect_<N>_opt1 (relax)
        defect_<N>_stat (static)
    {end}

    {begin}
    defect_<B>_stat (static_to_neb), defect_<E>_stat (static_to_neb)
        neb_<B-E>_opt1 (neb)
    {end}

**The <Q> tag**
The <Q> tag will correspond to charges given in :doc:`3_1_5_defects`.

* Charges are given as 
    
    * q=p0 for no charge
    * q=nX for negative charge X (addition of electrons)
    * q=pX for positive charge X (removal of electrons)

**{begin} and {end}**

In the recipe, {begin} and {end} will loop over, match up, and fill in scaling labels <S>, defect labels <N, B, and E>, NEB labels <B-E>, charges <Q>, and phonons <P>

* Only charges in the charge range of both the <B> and <E> defect parents of an NEB will produce an charged NEB.

* Use a new {begin} and {end} when you have a new tree branch or unindentation in the recipe that switches between <N> and <B> or <E>

* Note that defect endpoints need to be the parents of all NEB optimizations and NEB static calculations. Therefore, the endpoint-neb parent-child block may look like the following::

    {begin}
    defect_<B>_stat (static_to_neb), defect_<E>_stat (static_to_neb)
        neb_<B-E>_opt1 (neb)
            neb_<B-E>_opt2 (neb)
                neb_<B-E>_stat (neb_static)
        neb_<B-E>_opt2 (neb)
        neb_<B-E>_stat (neb_static)
    {end}


Full example::

    $recipe
    perfect_opt1 (lowmesh)
        perfect_opt2
            perfect_stat (static)
            {begin}
            inducescaling_<S>
                inducedefect_<S>_<N> (inducedefect)
                    defect_<S>_<N>_<Q>_opt1 (lowmesh_defect)
                        defect_<S>_<N>_<Q>_opt2 (defect_relax)
                            defect_<S>_<N>_<Q>_stat (static)
            {end}
    {begin}
    defect_<S>_<N>_<Q>_stat (static)
        phonon_<S>_<N>_<Q>_<P> (phonon)
    {end}
    {begin}
    defect_<S>_<B>_<Q>_stat (static_to_neb), defect_<S>_<E>_<Q>_stat (static_to_neb)
        neb_<S>_<B-E>_<Q>_opt1 (neb_to_neb)
            neb_<S>_<B-E>_<Q>_opt2 (neb_to_nebstat)
                neb_<S>_<B-E>_<Q>_stat (nebstat_to_phonon)
        neb_<S>_<B-E>_<Q>_opt2 (neb_to_nebstat)
        neb_<S>_<B-E>_<Q>_stat (nebstat_to_phonon)
    {end}
    {begin}
    neb_<S>_<B-E>_<Q>_stat (nebstat_to_phonon)
        phonon_<S>_<B-E>_<Q>_<P> (phonon)
    {end}
    $end

.. raw:: html

    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-54660326-1', 'auto');
      ga('send', 'pageview');

    </script>

