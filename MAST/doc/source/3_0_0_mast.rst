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
