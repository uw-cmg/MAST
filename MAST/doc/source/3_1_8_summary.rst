###################################
The Summary section
###################################

The ``$summary`` section of the input file will cause a ``SUMMARY.txt`` file to be printed into the recipe directory, once the recipe is complete.

Each line in the summary section follows the format::

    <ingredient name search string> <summary information>

*  <search string> is a search string for matching ingredient names.

*  <summary information> refers to a python file in ``<MAST installation directory>/summary`` which is supposed to extract information from a given ingredient directory.

*  For example, the following section would extract energies from each ingredient matching ``vac`` in its name.

    $summary
    vac energy
    $end

