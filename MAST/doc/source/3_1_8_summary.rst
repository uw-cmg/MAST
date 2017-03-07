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


.. raw:: html

    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-54660326-1', 'auto');
      ga('send', 'pageview');

    </script>

