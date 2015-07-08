########################
Trying out MAST
########################
#.  Go to ``$HOME/MAST/examples`` (or ``$WORK/MAST/examples`` or a similar folder, if you moved the ``$HOME/MAST`` folder from its default location.)
#.  Select one of the examples. The fastest one is ``simple_optimization.inp``
#.  Copy that file::

        cp simple_optimization.inp test.inp

#.  Modify the test.inp file with the correct ``mast_exec``, ``mast_ppn``, ``mast_queue``, ``mast_walltime``, and other settings described in :doc:`Input File<3_0_inputfile>`

#.  Try to parse the input file, entering the following command as one line::

        nice -n 19 mast -i test.inp 

    *  The ``nice -n 19`` keeps this command low priority, since it is being run on the headnode (but it is not too intensive).
    *  The ``-i`` signals to MAST that it is processing an input file.
#. Your ``$MAST_SCRATCH`` directory should now have a recipe directory in it.

    * The recipe directory will have a name corresponding to the elements and the input file, and ending with a timestamp of YYYYMMDD"T"hhmmss. 
    * The recipe directory will contain several subfolders, which are ingredient directories.
#. Go to that recipe directory.

    *  To see the input options:

        *  ``cat input.inp`` (should be identical to test.inp since no looping was used)
        
            *  Note that you can use other viewing commands, not just ``cat``, but be careful not to edit any of these files.

        *  ``cat archive_input_options.txt`` (should show Al instead of element X1)
    *  To see information about the ingredient relationships MAST detected from the recipe template:

        *  ``cat archive_recipe_plan.txt``
        
        *  Look at the ``$personal_recipe`` section in the ``input.inp`` file
    
    *  To see ingredient statuses at a glance:

        *  ``cat status.txt``

#.  Run mast once: ``nice -n 19 mast``

#.  You should see a "mastmon" job appear on the queue specified in ``$MAST_CONTROL/mastmon_submit.sh``

#.  MAST should have detected that the first ingredient was ready to run, so when that process disappears, run mast again: ``nice -n 19 mast``

#.  Now you should see ``perfect_opt1`` appear on the queue.

#. ``status.txt`` in the recipe directory in ``$MAST_SCRATCH`` should show that ``perfect_opt1`` has a status of "Proceed to Queue", or "P".

#.  When the queued ``perfect_opt1`` job starts running, you should be able to see output files inside ``$MAST_SCRATCH/<recipe directory>/perfect_opt1``

#.  If you forgot some step above, or you encounter some errors, remove the recipe folder from ``$MAST_SCRATCH`` and start again from the beginning of this section.

#.  The ``$MAST_CONTROL`` folder gives you error messages and other information. See :doc:`Running MAST <5_0_runningmast>` for tips.

----------------------
Unit testing
----------------------

Unit tests are available through the MAST tar.gz file. See :doc:`12_0_programming`. (Unit tests are not installed by default using pip.)

To check the validity of the MAST source code, navigate to ``<MAST installation directory>/MAST/test`` and run the unit tests with::

    nosetests --exe

Some tests may have been designated to be skipped.
Errors should be reported to the MAST development team as an issue on the github site (see :doc:`Programming for MAST <12_0_programming>`).
