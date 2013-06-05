=================
Repetitions
=================
For you may want to do several repetitions of a recipe with minor changes.
Keywords ``indeploop`` and ``deploop`` in the input file will make multiple recipes.

In the case of ``indeploop``, the recipes will be independent of each other::
    
    indeploop sys Al, Cr, Fe

In the case of ``deploop``, the recipes will be dependent on each other in sequence. They will require a python class for evaluating the dependency (how each recipe starts and ends). Override methods in the recipebatch.py file (RecipeBatch class) in order to customize behavior::

    deploop recipebatch.py mast_kpoints 2x2x2 M, 4x4x4 M, 6x6x6 M

Only the following keywords are supported for looping:

* indeploop
    * sys

* deploop
    * mast_kpoints

