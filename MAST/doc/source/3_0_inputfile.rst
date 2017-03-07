###################################
Input File
###################################

When you use the command ``mast -i <inputfile>.inp``, MAST does the following:

    * Reads the input file
    
    * Creates a recipe directory in ``$MAST_SCRATCH``
    
    * Creates the ingredient directories under that recipe directory
    
    * Creates all the necessary metadata.txt files for that recipe and its ingredients.

MAST will then copy the input file into that recipe directory, as ``input.inp``.

MAST will refer to this recipe-local ``input.inp`` file for all subsequent contact with the recipe.

*************************************
General structure of the input file
*************************************

The input file has many sections. Sections are denoted by ``$<section name>`` and ``$end``::

    $section
    
    section_text
    section_keyword keyword_value

    $end

Within each section there may also be subsections, with keywords and values.
Subsections are denoted by ``begin <subsection name>`` and ``end``.::

    $section
    
    section_text
    section_keyword keyword_value

    begin subsection
    subsection_text
    subsection_keyword subsection_keyword_value
    end

    $end

* Comments in the input file are allowed only as separate lines, starting with the # sign. 
* A comment may not be appended to a line.

*****************************
Sections of the input file
*****************************
See :doc:`3_1_inputsections`

******************************
Looping in the input file
******************************
If special looping tags are present in the input file, MAST can read in a single input file and create several permutated recipes in ``$MAST_SCRATCH``.

The looping tag ``indeploop`` may be used to create combinatorial permutations.

* ``indeploop`` may be used once at the beginning of a line (that is not a section or subsection header or "end" line).

* ``indeploop`` may be used multiple times in an input file.

When ``indeploop`` is present at the beginning of the line, input file permutations will be created depending on the values in parentheses. ::

    indeploop keyword1 (k1value1,k1value2)

The previous line would create two input files and corresponding recipes.
On the line that used to have ``indeploop`` on it, one input file would have::

    keyword1 k1value1

The other input file would have::

    keyword1 l1value2

If indeploop tags are present multiple times in the recipe, input files are created combinatorially::

    indeploop keyword1 (k1value1,k1value2)
    indeploop keyword2 (k2value1,k2value2)

The previous two lines in an input file would create four input files and corresponding recipes.
One input file would have::

    keyword1 k1value1
    keyword2 k2value1

Another would have::

    keyword1 k1value1
    keyword2 k2value2

A third would have::

    keyword1 k1value2
    kewyord2 k2value1

A fourth would have::

    keyword1 k1value2
    keyword2 k2value2

Sometimes, instead of combinatorial looping, some loops are meant to go together. In this case, the ``pegloop1`` and ``pegloop2`` tags may be used.

* There are only two pegged looping tags allowed, ``pegloop1`` and ``pegloop2``.
* Each tag may be used only once on a line.
* Each tag may be used on multiple lines.

Every line that starts with ``pegloop1`` (the same will apply for ``pegloop2``) will loop over keyword values, much like ``indeploop``. However, the point of the pegged loops is to have two or more keywords loop together.

For example::

    pegloop1 keyword1 (k1value1,k1value2)
    pegloop1 keyword2 (k2value1,k2value2)

Using the ``pegloop1`` tag, the lines above would not produce four input files and corresponding recipes, as they would when using the ``indeploop`` tag. Instead, they would produce only two input files and corresponding recipes.

One input file would have::

    keyword1 k1value1
    keyword2 k2value1

The other input file would have::

    keyword1 k1value2
    keyword2 k2value2

The number of items in parentheses should be equal for all instances of the ``pegloop1`` (or, separately, the ``pegloop2``) tag.

``pegloop1``, ``pegloop2``, and all instances of ``indeploop`` will work combinatorially with each other.

Complex example (for looping only - many other necessary lines in the input file are skipped)::    

    $mast
    pegloop1 system_name (strain1,strain2,strain3)
    $end

    $structure
    begin lattice
    pegloop1 (3,4,5) 0 0
    pegloop1 0 (3,4,5) 0
    pegloop1 0 0 (3,4,5)
    end

    begin elementmap
    pegloop2 X1 (Cr,Mn)
    end
    $end

    $ingredients
    begin ingredients_global
    indeploop mast_xc (pw91,pbe)
    LDAUJ 1
    pegloop2 LDAUU (4.5,5)
    end
    $end

The above example would create 3*2*2 = 12 input files and corresponding recipes. The input file for the one of the recipes would look like::

    $mast
    system_name strain2
    $end

    $structure
    begin lattice
    4 0 0
    0 4 0
    0 0 4
    end

    begin elementmap
    X1 Mn
    end
    $end

    $ingredients
    begin ingredients_global
    indeploop mast_xc pbe
    LDAUJ 1
    LDAUU 5
    end
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

