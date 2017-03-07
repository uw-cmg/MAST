###############################
Old Installation Instructions
###############################
 
=================================
Do pre-installation steps
=================================

-----------------------------
Locate your user profile
-----------------------------

Your user profile will set up environment variables like ``$PATH`` when you log in.

This installation will ask you to modify your user profile several times.

If you are comfortable modifying your user profile, please skip to :ref:`old-use-cluster-correctly`.

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

.. _old-use-cluster-correctly:

------------------------------
Use your cluster correctly
------------------------------
For this installation, please follow the correct procedures in order to avoid excessive headnode use on your cluster.

*  For example, you may want to preface every command with ``nice -n 19`` in order to reduce headnode load. 

*  Or, your cluster may have a dedicated compile node, or it may support interactive queue submission.

Please check with your cluster administrator if you are unsure of the correct procedures to user.

================================
Verify your Python version
================================

Check your version of python: ``python --version``

*  If your version of python is a 2.7 version (e.g. 2.7.3), skip to :ref:`verify-numpy-scipy`.

*  Otherwise, go to :ref:`locate-other-python`.

.. _locate-other-python:

-------------------------------------------------------------------
Locate and use an available but non-default 2.7 version of python
-------------------------------------------------------------------

For clusters using the "module" system, like Stampede or DLX, check which modules are available using ``module avail python`` or ``module avail``

For clusters not using the module system, you may need to look in ``//share/apps`` or in a similar shared directory, or ask your system administrator.

*  If you cannot find an existing version of python 2.7, skip to :ref:`install-local-python`.

*  If you did find an existing version of python 2.7, make sure that it is defaulted to be used first.

    #.  Add the appropriate line to your user profile. 

        *  If the version found was a module, then:
        
            *  Add a line like the following::

                module load python

            *  The actual wording may depend on the module name. Here, we are assuming that "python" is the module name.

        *  If the version found was located in an explicit directory and not found through the module system:
        
            *  Add a line like the following::

                export PATH=<path_to_python2.7>:$PATH

            *  For example::

                export PATH=//share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin:$PATH

    #.  Then, log out and log back in. 

    #.  Type ``which python`` or ``python --version`` to make sure your default version is now the correct version.

If you already use python for something else and shifting python versions will interfere with other programs, for example, you routinely use Python 2.4.3 instead and your other programs break if called from python 2.7.3, please contact your system administrator or the MAST development team.

Otherwise, go on to :ref:`verify-numpy-scipy`.

.. _verify-numpy-scipy:

-----------------------------------------
Verify that python has numpy and scipy
-----------------------------------------

Check to see if your python version has numpy and scipy::

    python

And then, from the python prompt::

    import numpy
    import scipy

If you receive an ImportError, then you must install a local version of python which has numpy and scipy. Go to :ref:`install-local-python`.

.. _install-local-python:

------------------------------------------------------------
Install a local version of python with numpy and scipy
------------------------------------------------------------

The EPD/Canopy version is preferred because it includes numpy and scipy already. Download this version from `EPD Free Canopy <https://www.enthought.com/downloads/>`_

*  Run the setup script. (e.g. ``bash ./canopy-1.0.3-rh5-64.sh``)

*  Follow the prompts and specify a local installation (use spacebar to scroll through the license file).

Add a line to your user profile to make this python installation your default python, for example::

    export PATH=//home/<username>/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin:$PATH

*  Do not just use the ``Canopy/bin`` directory, as python modules will not load properly

*  Log out and log back in.

Check your version of python: ``python --version``

The version given must be the correct version.

Check that numpy and scipy are installed, which they should be::

    python

And then at the python prompt::

    import numpy
    import scipy

==============================================
Install dependencies
==============================================
MAST requires pymatgen and custodian, each of which has several dependencies, which also have their own dependencies.

---------------------------------
pip note for the python-savvy
---------------------------------
If you have ``pip``, it is possible but sometimes unusually complicated to use pip to install MAST and its dependencies. 

*  If the ``pip`` command does not exist (``which pip`` does not return anything), go on to :ref:`manual-installation`.

If you have the ``pip`` command, it may be worth trying the following::

    pip install pymatgen==2.8.7 --user
    pip install custodian==0.7.5 --user
    pip install MAST --user

*  If this series of commands actually worked without errors, then do a quick installation of ASE following the instructions on the `ASE website <https://wiki.fysik.dtu.dk/ase/download.html>`_ and then skip to :ref:`add-local-bin`. 

*  If you have never used pip before, and using pip created a ``$HOME/.local`` folder for you for the first time, and you encounter errors, delete the ``$HOME/.local`` folder and go on to :ref:`manual-installation`.

*  If you encountered errors and your ``$HOME/.local`` folder already existed, carefully remove the most recent package folders under ``$HOME/.local/lib/python2.7/site-packages`` and go on to :ref:`manual-installation`.

.. _manual-installation:

----------------------------------------
Install dependencies manually
----------------------------------------

Download ``tar.gz`` files for the following dependencies from the `Python Package Index <https://pypi.python.org>`_

*  The versions listed are known to be compatible with MAST and with each other.

*  Using other version numbers may require adjustments to the entire list.

    *  In this case, look at ``install_requires`` inside the ``setup.py`` file to see which version numbers may be required.

Dependency list::

    PyCifRW-3.6.2.tar.gz
    pybtex-0.18.tar.gz
    pyhull-1.4.5.tar.gz
    monty-0.3.4.tar.gz
    PyYAML-3.11.tar.gz
    requests-2.3.0.tar.gz
    pymatgen-2.8.7.tar.gz
    custodian-0.7.5.tar.gz

Also get::

    python-ase-3.8.1.3440.tar.gz
    
from the `ASE website <https://wiki.fysik.dtu.dk/ase/download.html>`_

Upload each of these .tar.gz files onto your cluster.
Uncompress and untar each of these files (``tar -xzvf <tar.gz filename>``, for example, ``tar -xzvf PyCifRW-3.6.2.tar.gz``).

Following the order listed above, go to the untarred directory for each package and run the setup script as follows::

    tar -xzvf PyCifRW-3.6.2.tar.gz
    cd PyCifRW-3.6.2
    python setup.py install (--user, depending on the notes below)

And so on for all the packages.

If you are using a system-wide python, like from the module system or in a shared directory, then you need the ``--user`` tag, and will use the command::
    
    python setup.py install --user
    
In this case, the modules will end up in a folder like ``//home/<username>/.local/lib/python2.7/site-packages``.

If you are using your own locally-installed python, you can just use::

    python setup.py install
    
In this case, the modules will end up in your python installation directory, for example, ``//home/<username>/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/lib/python2.7/site-packages``.

If pymatgen cannot be installed because gcc cannot be found in order to compile spglib, then please see your system administrator.


If you will be using the mast_DFE_tool, then you will also need::

    pip install numpy --upgrade
    pip install pandas
    
.. _add-local-bin:

------------------------------------------------
Add the .local/bin directory, if necessary
------------------------------------------------

If you have a ``$HOME/.local/bin`` directory from a ``--user`` installation from any of the previous steps, add this directory to your ``$PATH`` environment variable by adding a line to your user profile, for example::
    
    export PATH=$HOME/.local/bin:$PATH

(This line can go either before or after any other ``export PATH`` lines you might have in your user profile.)

Then log out and log back in.

If you were using your own locally-installed python, then you would have already added the correct bin directory to your user profile in the :ref:`install-local-python` step. 


.. _old-vasp-psp-dir:

======================================
Set up the pymatgen VASP_PSP_DIR
======================================
This step is necessary if you are running VASP with MAST. If you are not running VASP with MAST, skip to :ref:`install-mast`.

--------------------------------------
Set up the pseudopotential folders
--------------------------------------

Locate the VASP pseudopotentials. If you cannot locate the VASP pseudopotentials, contact your system administrator or another person who uses VASP on the cluster.

``which potcar_setup.py`` should return the pymatgen utility for setting up your pseudopotential directories in the way that pymatgen requires.
If this command does not return a file location, then probably ``$HOME/.local/bin`` or ``<python installation directory>/bin`` is missing from your ``$PATH`` environment variable. See :ref:`add-local-bin`.

Run ``potcar_setup.py``::

    potcar_setup.py

*  The first directory address that you give to the utility is the directory that contains a few subdirectories, for example: potpaw_GGA, potpaw_LDA.52, potpaw_PBE.52, potUSPP_LDA, potpaw_LDA, potpaw_PBE, potUSPP_GGA. 

    *  These subdirectories themselves contain many sub-subdirectories with element names like Ac, Ac_s, Zr_sv, etc.

*  The second directory address that you give should be a new directory that you create.

Once the new pymatgen-structured folders have been created, rename the GGA PBE folder to ``POT_GGA_PAW_PBE``.

Later on, ingredients with a value of ``pbe`` for the ingredient keyword ``mast_xc`` will draw pseudopotentials out of this folder (see :doc:`3_0_inputfile`). 

Rename the GGA PW91 folder to ``POT_GGA_PAW_PW91``. Ingredients with a value of ``pw91`` for the ingredient keyword ``mast_xc`` will draw pseudopotentials out of this folder.

Example of running the python setup tool::
        
    Please enter full path where the POT_GGA_PAW_PBE, etc. 
    subdirs are present. 
    If you obtained the PSPs directly from VASP, this should 
    typically be the directory that you untar the files to : 
    
    //share/apps/vasp_pseudopotentials/paw
    
    Please enter the fullpath of the where you want to create 
    your pymatgen resources directory:

    //home/<username>/.local/vasp_pps

Rename the folders under ``//home/<username>/.local/vasp_pps``::

    mv //home/<username>/.local/vasp_pps/<pbe_name> //home/<username>/.local/vasp_pps/POT_GGA_PAW_PBE

    mv //home/<username>/.local/vasp_pps/<pw91_name> //home/<username>/.local/vasp_pps/POT_GGA_PAW_PW91

For assistance with potcar_setup.py, please contact the
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


.. _install-mast:

===============================
Install MAST
===============================
(If you successfully used ``pip`` to install MAST, go to :ref:`old-mast-setup`.)

* Get the `latest MAST package from the Python Package Index <https://pypi.python.org/pypi/MAST>`_

* Extract the package using ``tar -xzvf MAST-<version number>.tar.gz``

* Change into the package directory and run ``python setup.py install`` or ``python setup.py install --user`` as you did with the other packages in :ref:`manual-installation`.

You should be prompted to set the MAST environment variables, which is covered in :ref:`old-mast-setup`.

.. _old-mast-setup:

-----------------------------------------------------
Set the MAST environment variables
-----------------------------------------------------

The MAST setup.py script should have set up a ``MAST`` directory in your home directory, that is, ``//home/<username>/MAST``.

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

    *  Go to :ref:`old-additional-setup`.

*  If your platform is not matched exactly, or you would choose one of the generic choices:

    *  Set the three other environment variables (MAST_SCRATCH, MAST_ARCHIVE, and MAST_CONTROL) in your user profile.
    
    *  Log out and log back in.
    
    *  Go to :ref:`old-make-custom-platform`.

.. _old-make-custom-platform:

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

.. _old-additional-setup:

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

