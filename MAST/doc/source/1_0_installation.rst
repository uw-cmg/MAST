#############
Installation
#############
 
=================================
Do pre-installation steps
=================================

-----------------------------
Locate your user profile
-----------------------------

Your user profile will set up environment variables like ``$PATH`` when you log in.

This installation will ask you to modify your user profile several times.

If you are comfortable modifying your user profile, please skip to :ref:`use-cluster-correctly`.

For others: 

*  Your user profile is probably located in your home directory as ``//home/<username>/<user profile name>``, for example, ``//home/<username>/.bashrc``

*  Common user profile names are ``.bashrc``, ``.bash_profile``, ``.profile``, and ``.profile_user``

    *  These names usually start with a dot. 
    
    *  You may need to use the command ``ls -a`` to see these "hidden" files.

    *  Sometimes you may need to create your own user profile file. 
    
        *  For example, you may have a ``.profile`` file listed, but when you look at it, it tells you to create and modify a ``.profile_user`` file.

*  **After you save your changes to the user profile, you need to log out and then log back in, in order to see the changes take effect.**

    *  Alternately, you can ``source <user profile name>``, but occasionally this command will produce complications, for example, in path order.

If you cannot locate your user profile, please contact your system administrator.

.. _use-cluster-correctly:

------------------------------
Use your cluster correctly
------------------------------
For this installation, please follow the correct procedures in order to avoid excessive headnode use on your cluster.

*  For example, you may want to preface every command with ``nice -n 19`` in order to reduce headnode load. 

*  Or, your cluster may have a dedicated compile node, or it may support interactive queue submission.

Please check with your cluster administrator if you are unsure of the correct procedures.

================================================
Install local python version and dependencies
================================================

MAST requires the following dependencies, some of which have additional dependencies:

*  numpy
*  scipy
*  matplotlib
*  pymatgen
*  custodian
*  pandas
*  ase

Currently, the Anaconda python package is preferred because it:

*  automatically installs to a user directory and does not need root privileges
*  does not have the environment complications that the previously-recommended Enthought Canopy python has
*  can easily install numpy, scipy, and matplotlib.


-------------------------------
Install Anaconda python
-------------------------------

Download the free installer from `Anaconda <https://store.continuum.io/cshop/anaconda/>`_. Use a Python 2.X version.

*  Run the setup script. (e.g. ``bash ./Anaconda-<version>.sh``)

*  Follow the prompts and specify a local installation (use spacebar to scroll through the license file).

When prompted, agree for the installer to add a line to your user profile to make this python installation your default python, for example::

    export PATH=//home/<username>/anaconda/bin:$PATH

*  Log out and log back in. Verify that your python version is now the anaconda python version.::

    which python

Use the conda command to install the following packages::

    conda install numpy scipy matplotlib

If you are developing for MAST, install the following packages::

    conda install numpy scipy matplotlib nose sphinx

--------------------------------------------
Install the MAST package
--------------------------------------------

Install the MAST package from the python package index::

    pip install MAST

This command should install MAST, pymatgen, custodian, and pandas.

If pip cannot find anything (SSL certificate fail), try::

    pip install MAST --cert <path to conda>/ssl/cacert.pem

If pymatgen cannot be installed because gcc cannot be found in order to compile spglib, then please see your system administrator.

If you need additional standalone tools (see :doc:`8_0_standalonetools`) or unit tests, then also get and unzip the MAST tar.gz file (see :doc:`12_0_programming`).


.. _vasp-psp-dir:

======================================
Set up the pymatgen VASP_PSP_DIR
======================================
This step is necessary if you are running VASP with MAST. If you are not running VASP with MAST, skip to :ref:`install-mast`.

--------------------------------------
Set up the pseudopotential folders
--------------------------------------

Locate the VASP pseudopotentials. If you cannot locate the VASP pseudopotentials, contact your system administrator or another person who uses VASP on the cluster.

Follow the POTCAR Setup instructions from Pymatgen's installation instructions at `http://pymatgen.org <http://pymatgen.org>`_

Later on, ingredients with a value of ``pbe`` for the ingredient keyword ``mast_xc`` will draw pseudopotentials out of the ``POT_GGA_PAW_PBE`` folder (see :doc:`3_0_inputfile`). 

Ingredients with a value of ``pw91`` for the ingredient keyword ``mast_xc`` will draw pseudopotentials out of the ``POT_GGA_PAW_PW91`` folder, and so on. 

Use pymatgen-recognized values for ``mast_xc`` keyword values.

For assistance with POTCAR setup, please contact the
`Pymatgen support group <https://groups.google.com/forum/#!forum/pymatgen>`_

---------------------------------------------
Add the VASP_PSP_DIR to your user profile
---------------------------------------------
Add a line to your user profile exporting the environment variable ``$VASP_PSP_DIR`` to the new pseudopotential directory created above.

For example::

    export VASP_PSP_DIR=//home/<username>/.local/vasp_pps

Log out and log back in.

Test the change::
    
    cd $VASP_PSP_DIR

*  Make sure you are getting to the right directory, which has the ``POT_GGA_PAW_PBE`` etc. folders inside it.



.. _mast-setup:

======================================
Set the MAST environment variables
======================================

The pip installation of MAST should have set up a ``MAST`` directory in your home directory, that is, ``//home/<username>/MAST``.

**With conda installation** it appears that this folder is buried under your anaconda directory, like ``./anaconda/lib/python2.7/site-packages/home/<username>/MAST``. Go ahead and move that MAST folder into your home directory::

    mv ~/anaconda/lib/python2.7/site-packages/home/<username>/MAST ~/.

*  This directory is primarily for storing calculations, and should not be confused with the python module directory, which is where the actual MAST python code resides.

Inside ``$HOME/MAST`` there should be:

#.  A ``SCRATCH`` folder:

    *  Each time an input file is given to MAST, MAST will create a recipe directory inside this folder. 
    
    *  Each recipe directory will itself contain ingredient, or calculation, directories. Calculations will be submitted to the queue from inside these ingredient directories. 

    *  Multiple recipes may reside in ``SCRATCH`` at the same time, and MAST will evaluate them alphabetically.

#.  An ``ARCHIVE`` folder: 

    *  When a recipe directory is complete, MAST will move it from ``SCRATCH`` to ``ARCHIVE``.

#.  A ``CONTROL`` folder: 

    *  MAST requires some control files in order to run. It also does some higher-level logging, and stores that output here.

*  On some clusters, like Stampede, the home directory is not where you actually want to store calculations. Instead, there may be a separate "work" or "scratch" directory. In this case, move the entire ``$HOME/MAST`` directory into the work or scratch directory, for example::

    mv $HOME/MAST $WORK/.

In this case, the environment variables below should therefore say ``$WORK`` instead of ``$HOME``.

*  You can also move the MAST directory anywhere else, as long as you set the environment variables correcty.

Copy and paste the environment variables into your user profile, setting the paths correctly if you have moved the ``$HOME/MAST`` directory::

    export MAST_SCRATCH=$HOME/MAST/SCRATCH
    export MAST_ARCHIVE=$HOME/MAST/ARCHIVE
    export MAST_CONTROL=$HOME/MAST/CONTROL
    export MAST_PLATFORM=<platform_name>

For platform_name, choose from one of the following::
    
    aci
    bardeen
    dlx
    korczak
    no_queue_system
    pbs_generic
    sge_generic
    slurm_generic
    stampede
    turnbull

For example::

    export MAST_PLATFORM=stampede

*  If your platform is available by name (not _generic), then:

    *  Add the four environment variable lines to your user profile as above.
    
    *  Log out and log back in.

    *  Go to :ref:`additional-setup`.

*  If your platform is not matched exactly, or you would choose one of the generic choices:

    *  Set the three other environment variables (MAST_SCRATCH, MAST_ARCHIVE, and MAST_CONTROL) in your user profile.
    
    *  Log out and log back in.
    
    *  Go to :ref:`make-custom-platform`.

.. _make-custom-platform:

---------------------------------------
Make a custom platform, if necessary
---------------------------------------

Run the following command. It should produce some errors, but ignore those and just see where MAST is installed::

    mast -i none

For example, output may be::

    ------------------------------------------------------
    Welcome to the MAterials Simulation Toolkit (MAST)
    Version: 1.1.5
    Installed in: .local/lib/python2.7/site-packages/MAST
    ------------------------------------------------------

and then some errors.

Go to the "installed in" directory, and then::

    cd submit/platforms

Identify the closest-matching directory to your actual platform (for example, if you have an SGE platform, this directory would be sge_generic)

Copy this directory into a new directory inside the ``platforms`` folder, for example::

    cp -r sge_generic my_custom_sge

Then, inside your new folder, like ``my_custom_sge``, modify each of the following files as necessary for your platform::

    submit_template.sh
    mastmon_submit.sh
    queue_commands.py

Explanations for each file are given in the following sections. Modify and test each file in your new custom platform folder.

Then, in your user profile, use your new custom folder for the platform name of ``$MAST_PLATFORM``::

    export MAST_PLATFORM=my_custom_sge

Log out and log back in.

^^^^^^^^^^^^^^^^^^^^^^^^
submit_template.sh
^^^^^^^^^^^^^^^^^^^^^^^^

``submit_template.sh`` is the generic submission template from which ingredient submission templates will be created.

*  MAST will replace anything inside question marks, for example ``?mast_ppn?`` with the value of the appropriate keyword.

The following keywords may be used; see :doc:`Input File <3_0_inputfile>` for more information on each keyword.

* mast_processors
* mast_ppn
* mast_nodes
* mast_queue
* mast_exec
* mast_walltime
* mast_memory
* mast_name (the ingredient name)

Examine the template carefully, as an error here will prevent your ingredients from running successfully on the queue.

*  The provided template should be a good match for its platform.

    *  Otherwise, you can take one of your normal submission templates and substitute in ``?mast_xxx?`` fields where appropriate.

*  Or, vice versa, you can take the provided template, replace the ``?mast_xxx?`` fields with some reasonable values, and see if filled-in submission template will run a job if submitted normally using ``qsub``, ``sbatch``, etc.

^^^^^^^^^^^^^^^^^^^^^^^^
mastmon_submit.sh
^^^^^^^^^^^^^^^^^^^^^^^^

``mastmon_submit.sh`` is the submission template that will submit the MAST Monitor to the queue every time ``mast`` is called.

The MAST Monitor will check the completion status of every recipe and ingredient in the ``$MAST_SCRATCH`` folder.

*  If you have a recipe you would like to skip temporarily, manually put a file named ``MAST_SKIP`` inside that recipe's folder in ``$MAST_SCRATCH``. ``MAST_SKIP`` can be an empty file, or it can contain notes; MAST does not check its contents.

*  ``mastmon_submit.sh`` should be set to run on the shortest-wallclock, fastest-turnaround queue available, e.g. a serial queue

The ``mastmon_submit.sh`` script is copied into the ``$MAST_CONTROL`` directory every time you run ``mast``.

If you see that after you type ``mast``, no "mastmon" process appears on the queue, then test the submission script directly::

    cd $MAST_CONTROL
    qsub mastmon_submit.sh (or use sbatch for slurm, etc.)

*  Modify the ``$MAST_CONTROL/mastmon_submit.sh`` file until it the "mastmon" process successfully runs on the queue.

*  Copy your changes into the ``<MAST installation directory>/submit/platforms/<platform>/mastmon_submit.sh`` file so that your changes will be reflected the next time that you run MAST.

^^^^^^^^^^^^^^^^^^^^^^^^^^
queue_commands.py
^^^^^^^^^^^^^^^^^^^^^^^^^^

These queue commands will be used to submit ingredients to the queue and retrieve the job IDs and statuses of ingredients on the queue.

*  For a custom platform, modify the ``<MAST installation directory>/submit/platforms/<your custom platform>/queue_commands.py`` file.

*  Do not modify the ``<MAST installation directory/submit/queue_commands.py`` file.

Modify the following python functions as necessary:

*  ``queue_submission_command``: 

    *  This function should return the correct queue submission command, 
    
    *  For example, this function should return ``qsub`` on PBS/Torque, or ``sbatch`` on slurm.

*  ``extract_submitted_jobid``:

    *  This function should parse the job ID, given the text that returns to screen when you submit a job.
    
    *  For example, it should return ``456789`` as the jobid for the following job submission and resulting screen text::

        login2.mycluster$ sbatch submit.sh 
        -----------------------------------------------------------------
                  Welcome to the Supercomputer              
        -----------------------------------------------------------------
        --> Verifying valid submit host (login2)...OK
        --> Verifying valid jobname...OK
        --> Enforcing max jobs per user...OK
        --> Verifying job request is within current queue limits...OK
        Submitted batch job 456789

    *  On a different cluster, it would return ``456789`` as the jobid for the following submission and resulting screen text::

        [user1@mycluster test_job]$ qsub submit.sh
        456789.mycluster.abcd.univ.edu

*  ``queue_snap_command``:

    *  This function should show a summary of your current submitted jobs, which we call the ``queue_snapshot``.

    *  For example, the queue snapshot command should return something like the following (platform-dependent)::
    
        JOBID   PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
        456789      normal test1 user1 PD       0:00      4 (Resources)
        456788      normal test2 user1 PD       0:00      1 (Resources)
        456774      normal test3 user1  R    6:14:53      1 c123-124
        456775      normal test4 user1  R    6:15:34      1 c125-126

*  ``queue_status_from_text``:

    *  This function should return the status of a specific job, based on the job number.

    *  For example, job 456789 in the queue snapshot above, with status "PD" should correspond to a "Q" status (queued status) for MAST.

    *  Job 456775 in the queue snapshot above, with status "R", should correspond to an "R" status (running status) for MAST.

*  ``get_approx_job_error_file``:

    *  This function should return the name of the job error file.

    *  The name of this file will depend on what is specified in ``submit_template.sh`` and is usually something like ``slurm.<jobnumber>`` or ``<jobname>.e<jobnumber>``

.. _additional-setup:

================================
Additional setup
================================

You may need to do any or all of the following:

* Identify the correct ``mast_exec`` call for your system.

    *  For example, suppose you run VASP like this::

        //opt/mpiexec/bin/mpiexec //share/apps/bin/vasp5.2_par_opt1

    *  Then, in your input files, the ``mast_exec`` keyword would be specified like this::

        mast_exec //opt/mpiexec/bin/mpiexec //share/apps/bin/vasp5.2_par_opt1

*  Add additional lines to your user profile which allow you to run VASP, including any modules that need to be imported, additions to your library path, unlimiting the stack size, and so on.

*  Modify your text editor settings so that tabs become four spaces (or so that you have such an option readily available). This setting is very important to ensure that MAST can read the input file, especially the recipe section of the input file.

    *  If you use VIM (``vi``), add the following lines to your ``~/.vimrc`` file::
    
        " VIM settings for python in a group below:
        set tabstop=4
        set shiftwidth=4
        set smarttab
        set expandtab
        set softtabstop=4
        set autoindent

Once you have completed any additional setup and have identified what ``mast_exec`` should be, go to :doc:`17_0_testmast`.


.. raw:: html

    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-54660326-1', 'auto');
      ga('send', 'pageview');

    </script>

