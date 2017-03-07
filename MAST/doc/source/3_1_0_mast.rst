###################################
The MAST section
###################################
The ``$mast`` section contains this keyword:

**system_name**: Specify a single descriptive word here.

This keyword will become part of the recipe directory's name, and allow you to spot the recipe more easily in the ``$MAST_SCRATCH`` directory.

This keyword comes in handy with pegged looping, in order to help you identify loops.

    * Loops are otherwise differentiated by elements, if you were looping over elements, or simply by a 1-second timestamp difference.

Example::
    
    $mast
    system_name epitaxialstrain
    $end

Example for pegged loop::

    $mast
    pegloop1 system_name (strain1,strain2,strain3)
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

