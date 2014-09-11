#############
Installation
#############
 
*********************************
Installation
*********************************

===========================
Cluster Usage
===========================
Please follow the correct procedures to avoid excessive headnode use on your cluster.
For example, you may want to preface every command with ``nice -n 19`` to reduce headnode load. Or, your cluster may have a dedicated compile node or support interactive queue submission. Please check with your cluster administrator if you are not sure.

================================
Verify your Python version
================================
Check your version of python: ``python --version``

If your version of python is not a 2.7 version (e.g. 2.7.3), try to locate an existing version of python 2.7.

-------------------------------------------------------------------
Locating and using an available but non-default version of python
-------------------------------------------------------------------
For clusters using the "module" system, like Stampede or DLX, check which modules are available using ``module avail python``, or ``module avail``

For clusters not using the module system, you may need to look in ``//share/apps`` or a similar shared directory, or ask your system administrator.

If you cannot find an existing version of python 2.7, skip to :ref:`install-local-python`.

If you did find an existing version of python 2.7, now make sure that it is defaulted to be used first.
You may need to add a line to your user profile. 
Your user profile may be located in ``//home/username/.bashrc`` or a similar file.

If the version found was a module, add a line like::

    module load python

The actual wording may depend on the module name.

If the version found was found in an explicit directory and not through the module system, add a line like::

    export PATH=<path_to_python2.7>:$PATH

For example::

    export PATH=//share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin:$PATH

Then, log out and log back in. 

Type ``which python`` to make sure you have the right version, or ``python --version``.

If you already use python for something else and shifting python versions will interfere with other programs, for example, you routinely use Python 2.4.3 instead and your other programs break if called from python 2.7.3, please contact your system administrator or the MAST development team.

Now, check to see if your python version has numpy and scipy::

    python

And then, from the python prompt::

    import numpy
    import scipy

If you receive an ImportError, then you must install a local version of python which has numpy and scipy. Go to :ref:`install-local-python`.

.. _install-local-python:

------------------------------------------------------------
Installing a local version of python with numpy and scipy
------------------------------------------------------------

The EPD/Canopy version is preferred because it includes numpy and scipy already. Download this version from `EPD Free Canopy <https://www.enthought.com/downloads/>`_

*  Run the setup script. (e.g. ``bash ./canopy-1.0.3-rh5-64.sh``)
*  Follow the prompts and specify a local installation (use spacebar to scroll through the license file)

Add lines to your user profile to make this python installation your default python, for example::

    vi ~/.bashrc
    #EPD (Canopy) python
    export PATH=//home/<username>/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin:$PATH

*  Do not just use the ``Canopy/bin`` directory, as python modules will not load properly
*  Log out and log back in.

Check your version of python: ``python --version``

The version given must be the correct version.

Check that numpy and scipy are installed, which they should be::

    python
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

    pip install pymatgen --user
    pip install custodian --user
    pip install MAST --user

*  If this series of commands actually worked without errors, skip to :ref:`add-local-bin`. 
*  If you have never used pip before, and using pip created a ``$HOME/.local`` folder for you for the first time, and you encounter errors, delete the ``$HOME/.local`` folder and go on to :ref:`manual-installation`.
*  If you encountered errors and your ``$HOME/.local`` folder already existed, carefully remove the most recent package folders under ``$HOME/.local/lib/python2.7/site-packages`` and go on to manual installation.

.. _manual-installation:

----------------------------------------
Manual installation of dependencies
----------------------------------------

Download ``tar.gz`` files for the following dependencies from the `Python Package Index <https://pypi.python.org>`_

*  The versions listed are known to be compatible with MAST and with each other.

*  Using other version numbers may require adjustments to the entire list.
Look at ``install_requires`` inside the ``setup.py`` file to see which version numbers may be required.

Dependency list::

    PyCifRW-3.6.2.tar.gz
    pybtex-0.18.tar.gz
    pyhull-1.4.5.tar.gz
    monty-0.3.4.tar.gz
    PyYAML-3.11.tar.gz
    requests-2.3.0.tar.gz
    pymatgen-2.7.9.tar.gz
    custodian-0.7.5.tar.gz

Upload each of these .tar.gz files onto your cluster.
Uncompress and untar each of these files (``tar -xzvf <tar.gz filename>``, for example, ``tar -xzvf PyCifRW-3.6.2.tar.gz``).

In the order that they are given, go to the untarred directory for each package and run the setup script as follows::

    tar -xzvf PyCifRW-3.6.2.tar.gz
    cd PyCifRW-3.6.2
    python setup.py install (--user, depending on the notes below)

And so on for all the packages, in the order that they appear.

If you are using a system-wide python, like from the module system or in a shared directory, then you need the ``--user`` tag, and will use the command::
    
    python setup.py install --user
    
In this case, the modules will end up in a folder like ``//home/<username>/.local/lib/python2.7/site-packages``.

If you are using your own locally-installed python, you can just use::

    python setup.py install
    
In this case, the modules will end up in your python installation directory, for example, ``//home/<username>/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/lib/python2.7/site-packages``.

If pymatgen cannot be installed because gcc cannot be found in order to compile spglib, then please see your system administrator.

.. _add-local-bin:

------------------------------------------------
Adding the .local/bin directory, if necessary
------------------------------------------------

If you have a ``$HOME/.local/bin`` directory from a ``--user`` installation from any of the previous steps, add it to your PATH by modifying your user profile.
For example::
    
    vi ~/.bashrc
    export PATH=$HOME/.local/bin:$PATH

(This line can go either before or after any other ``export PATH`` lines you might have.)

Then log out and log back in.

If you were using your own locally-installed python, then you would have already added the correct bin directory to your user profile in the :ref:`install-local-python` step. 


.. _vasp-psp-dir:

======================================
Set up the pymatgen VASP_PSP_DIR
======================================
This step is necessary if you are running VASP with MAST. If you are not running VASP with MAST, you may skip this step.

Locate the VASP pseudopotentials. If you cannot locate the VASP pseudopotentials, contact your system administrator or another person who uses VASP on the cluster.

``which potcar_setup.py`` should return the pymatgen utility for setting up your pseudopotential directories in the way that pymatgen requires.

Usually, if pymatgen has been propery
set up the pymatgen VASP_PSP_DIR: 6.1 first find potcar_setup.py. usually if pymatgen has been properly installed, you can find potcar_setup.py by simply calling: which porcar_setup.py. However, in the rare case that porcar_setup.py is not in your path, you can find it under the installation path of pymatgen. Here’s how to find where pymatgen was installed: 1) find the paths for your site packages: >>> import site; site.getsitepackages() 2)check each of the site package path given above until finding the one containing a folder named pymatgen*.egg (e.g., pymatgen-2.10.2-py2.7-linux-x86_64.egg). 3). potcar_setup.py is in pymatgen*.egg/EGG-INFO/scripts/

6.2 run potcar_setup.py.
the “POT_GGA_PAW_PBE” directory address you give to potcar_setup.py is the directory that contains a few subdirectories, for example: potpaw_GGA potpaw_LDA.52 potpaw_PBE.52 potUSPP_LDA potpaw_LDA potpaw_PBE potUSPP_GGA. These subdirectory again contain many sub-subdirectory with element names like Ac, Ac_s, …, Zr_sv.

6.3 I found the following two sentences confusing:
• Rename the PBE folder POT_GGA_PAW_PBE to correspond to mast_xc pbe
• Rename the GGA folder POT_GGA_PAW_PW91 to correspond to mast_xc pw91
Now I got their meaning. I originally thought I need to rename “POT_GGA_PAW_PBE” to “mast_xc pbe”. So maybe we can add quotation marks and or change the font for “the PBE folder” and “POT_GGA_PAW_PBE”.
    
Run pymatgen's python setup tool. This tool should be located wherever pymatgen was installed, either ``~/.local/bin/potcar_setup.py`` if you installed it with ``--user``, or wherever python is, otherwise. ::

    python .local/bin/potcar_setup.py or python potcar_setup.py or simply potcar_setup.py
        
(Remember to use the correct version of python, determined in step 2, e.g. //share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python .local/bin/potcar_setup.py)

Take the paw directory if you are using PAW. Do not take the top directory, or the GGA/LDA/etc folders will overwrite.

Example of running the python setup tool::
        
    Please enter full path where the POT_GGA_PAW_PBE, etc. 
    subdirs are present. 
    If you obtained the PSPs directly from VASP, this should 
    typically be the directory that you untar the files to : 
    //share/apps/vasp_pseudopotentials/paw
    Please enter the fullpath of the where you want to create 
    your pymatgen resources directory:
    //home/<username>/.local/vasp_pps

Rename the folders under ``//home/<username>/.local/vasp_pps``:
    
*  Rename the PBE folder POT_GGA_PAW_PBE to correspond to mast_xc pbe
*  Rename the GGA folder POT_GGA_PAW_PW91 to correspond to mast_xc pw91

==============================================
Add the VASP_PSP_DIR to your user profile
==============================================
Add a line to your .bashrc file exporting the environment variable VASP_PSP_DIR to this VASP directory.

*  On bardeen, it should look something like::

    export VASP_PSP_DIR=//home/<username>/.local/vasp_pps

*  On DLX, use the directories already created::
    
    export VASP_PSP_DIR=//home/adozier/VASP/resources
    export VASP_PSP_DIR=<whichever path you used in the potcar_setup.py script>
*  Remember to save your .bashrc file. Test the change::
    
    source ~/.bashrc
    cd $VASP_PSP_DIR

*  Make sure you are getting to the right directory, which has POT_GGA_POW_PBE etc. folders inside it.

===============================
Install ASE
===============================

Obtain the latest source code from `<https://wiki.fysik.dtu.dk/ase/>`_

Unzip the tar.gz file to your home directory

Create the softlink as shown (use the version number you downloaded)::

    ln -s python-ase-3.8.0.3420 ase

DO NOT link to the ase folder within the unzipped tar.gz; only link to the top folder as shown above.

In your user profile, add the following line::

    export PYTHONPATH=$PYTHONPATH:~/ase

Log out and log back in.

===============================
Get MAST
===============================
* Get the latest MAST package from the Python package index::
    
    nice -n 19 pip install --upgrade --no-deps --user MAST

The no-dependencies tag is on because we are assuming pymatgen and custodian have been properly installed as above. It is recommended to install them separately.

Use the ``--user`` tag if you are not using the easy_install and pip from your own installation of python. Otherwise, you can omit this tag.
    
======================================
Set up the environment variables
======================================
The pip installation should set up a ``MAST`` directory in ``//home/username/MAST`` with several subdirectories.

The pip installation should then warn you with an ATTENTION flag of environment variables that must be set. 

You may copy and paste the environment variables from the terminal into your user profile. In the examples below, ``username`` should have been changed to your username.::
    
    export MAST_SCRATCH=//home/username/MAST/SCRATCH
    export MAST_ARCHIVE=//home/username/MAST/ARCHIVE
    export MAST_CONTROL=//home/username/MAST/CONTROL"
    export MAST_PLATFORM=platform_name

You will need to manually choose platform_name as one of the following::
    
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

    export MAST_PLATFORM=sge_generic

You must choose one of the platforms presented. Choose the best match. If your choice is not matched exactly, choose something anyway, complete the rest of this step, and go on to the following step.

Remember to log out and log back in after modifying your user profile.

-----------------------------------
Environment variable explanations
-----------------------------------
An explanation of each variable appears in the next section

MAST_SCRATCH: This variable may be set to any directory. MAST will look for recipes in this directory. ::
    
    export MAST_SCRATCH=//home/username/MAST/SCRATCH

MAST_ARCHIVE: This variable may be set to any directory. MAST will move completed recipes from ``$MAST_SCRATCH`` into this directory. ::
    
    export MAST_ARCHIVE=//home/username/MAST/ARCHIVE

MAST_CONTROL: This variable may be set to any directory. MAST monitor log files, MAST monitor error files, and other MAST monitor output will be written to this directory. ::
    
    export MAST_CONTROL=//home/username/MAST/CONTROL

MAST_CONTROL also has several subfolders. If you move your $MAST_CONTROL to a different path, please copy the subfolders with it.

MAST_PLATFORM: This variable switches among platforms. Note that it looks both in $MAST_CONTROL/platforms and in the platforms folder in your MAST installation directory (often in some path like //home/username/.local/lib/python2.7/site-packages/MAST or //share/apps/EPD.../lib/python2.7/site-packages/MAST). ::

    export MAST_PLATFORM=bardeen

VASP_PSP_DIR: This variable is necessary if VASP and VASP pseudopotential files are being used. See the documentation for the `Materials Project's <http://materialsproject.org>`_ `pymatgen <http://pymatgen.org>`_ code. The VASP_PSP_DIR should be set to a path which contains folder such as POT_GGA_PAW_PBE (for functional PBE, or mast_xc PBE in Ingredients) or POT_GGA_PAW_PW91 (for functional PW91). ::
    
    export VASP_PSP_DIR=//share/apps/MAST/vasp_pps

PATH: If you have created a local MAST installation using ``pip --install --no-deps --user``, then this variable should be appended with the ``//home/username/.local/bin`` directory so that the mast* executables may be found. ::
    
    export PATH=$PATH://home/username/.local/bin

Otherwise, if the mast executables are in ``//home/username/bin``, no such modification is needed.

=================================================
Modify submission details for your platform
=================================================
If your platform was not matched exactly, you or your system administrator should look where MAST was installed (e.g. often under some python folder, for example ``//share/apps/EPD...etc./lib/python-2.7/site-packages``, or, for a local installation, ``//home/username/.local/lib/python-2.7/site-packages``) and go to ``MAST/submit/platforms``.

Copy the closest-matching set of files into a new directory inside the ``platforms`` folder.
Then, modify each of the following files as necessary for your platform::

    submit_template.sh
    mastmon_submit.sh
    queue_commands.py

* Copy this new folder into your ``$MAST_CONTROL/platforms`` folder with the other platform folders.
* Edit ``$MAST_CONTROL/set_platform`` so that the word in it is the name of the new folder.
* Copy the new ``mastmon_submit.sh`` as ``$MAST_CONTROL/mastmon_submit.sh``

---------------------------------
mastmon_submit.sh
---------------------------------
This submission script is responsible for submitting to the ingredient- and recipe-checking script to the queue every time ``mast`` is called.

It should be set up to run on the shortest-wallclock, fastest-turnaround queue on your system (e.g. a serial queue, morganshort, etc.)

The script is copied into the $MAST_CONTROL directory by the ``initialize.py`` script and will be run from there.

Test mastmon_submit.sh by submitting it to the queue. A "mastmon" process should briefly appear on the queue. Continue to modify submit.sh until the "mastmon" process successfully runs on the queue.

Use commands similar to these (``sbatch`` instead of ``qsub`` for slurm)::

    cd $MAST_CONTROL
    qsub mastmon_submit.sh

------------------------------------------
submit_template.sh
------------------------------------------
This submission script template will be used to build submission scripts for the ingredients. Use ``?mast_keyword?`` to denote a place where the following MAST keywords (see :doc:`Input File <3_0_inputfile>` for more information on keywords) may be substituted in.

* mast_processors or a combination of mast_ppn and mast_nodes
* mast_queue
* mast_exec
* mast_walltime
* mast_memory
* the ingredient name

Examine the template carefully, as an error here will prevent your ingredients from running successfully on the queue.

-----------------------------
queue_commands.py
-----------------------------
These queue commands will be used to submit ingredients to the queue.

================================
Additional setup
================================
Figure out the correct mast_exec calls for your system, to be used in the :doc:`Input File<3_0_inputfile>`. Examples are below.

*  Bardeen: ``mast_exec //opt/mpiexec/bin/mpiexec //share/apps/bin/vasp5.2_par_opt1``  (or any of the other vasp executables) 
*  DLX: ``mast_exec //home/username/bin/vaspmpirun``, where vaspmpirun is the following script (indentations are all part of the previous line)::

    #!/bin/bash
    export PERL5LIB=/opt/moab/lib/perl5
    export MIC_LD_LIBRARY_PATH=/share/cluster/RHEL6.2/x86_64/
        apps/intel/ict/composer_xe_2013.0.079/compiler/lib/mic
    export LD_LIBRARY_PATH=/share/cluster/RHEL6.2/x86_64/apps/
        openmpi/1.6.2/lib:
        /share/cluster/RHEL6.2/x86_64/apps/intel/ict/
        composer_xe_2013.0.079/compiler/lib/intel64:
        /share/cluster/RHEL6.2/x86_64/apps/intel/ict/
        composer_xe_2013.0.079/mkl/lib/intel64
    export INTEL_MKL_LIBS=/share/cluster/RHEL6.2/x86_64/
        apps/intel/ict/composer_xe_2013.0.079/mkl/lib/intel64
    export QTLIB=/usr/lib64/qt-3.3/lib
    PATH=$PATH:$HOME/bin:$HOME/bin/convaspTest
    export PATH
    VaspPath=//home/adozier/VASP/vasp.5.2
    export OMP_NUM_THREADS=1
    ulimit -s unlimited
    ulimit -l unlimited
    #mpirun $VaspPath/vasp
    //share/cluster/RHEL6.2/x86_64/apps/openmpi/1.6.2/bin/
        mpirun $VaspPath/vasp

Modify ~/.bashrc if necessary
    
*  ACI/HPC, add: ``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH://opt/intel/lib/intel64``

To ensure recipes are created correctly, add python whitespace tab stops to your ~/.vimrc file::
    
    " VIM settings for python in a group below:
    set tabstop=4
    set shiftwidth=4
    set smarttab
    set expandtab
    set softtabstop=4
    set autoindent

Follow the testing instructions from :ref:`test-on-cluster`

.. _test-on-cluster:

*********************************
Test that MAST can run
*********************************
#.  Go to ``//home/username/MAST/examples``
#.  Select one of the examples. The fastest one is ``simple_optimization.inp``
#.  Copy that file::

        cp simple_optimization.inp test.inp

#.  Modify the test.inp file with the correct ``mast_exec``, ``mast_ppn``, ``mast_queue``, and other settings described in :doc:`Input File<3_0_inputfile>`

#.  Try to parse the input file, entering the following command as one line::

        nice -n 19 mast -i test.inp 

    *  The .nice -n 19. keeps this command low priority, since it is being run on the headnode (but it is not too intensive).
    *  The -i signals to MAST that it is processing an input file.
#. Your ``//home/username/MAST/SCRATCH`` directory should now have a recipe directory in it.

    * The recipe directory will have a name corresponding to the elements and the input file, and ending with a timestamp of YYYYMMDD"T"hhmmss. 
    * The recipe directory will contain several subfolders, which are ingredient directories.
#. Go to that recipe directory.

    *  To see the input options:

        *  ``cat input.inp`` (should be identical to test.inp since no looping was used)
        
            *  Note that you can use other viewing commands, not just .cat., but be careful not to edit any of these files.

        *  ``cat archive_input_options.txt`` (should show Al instead of element X1)
    *  To see information about the ingredient relationships MAST detected from the recipe template:

        *  ``cat personal_recipe.txt``
        *  ``cat archive_recipe_plan.txt``

    *  To see ingredient statuses at a glance:

        *  ``cat status.txt``

#.  Run mast once: ``nice -n 19 mast``
#.  You should see a `mastmon` job appear on the queue specified in $MAST_CONTROL/mastmon_submit.sh (which should be morganshort for bardeen).
#.  MAST should have detected that the first ingredient was ready to run, so when that process disappears, run mast again: ``nice -n 19 mast``
#.  Now you should see ``perfect_opt1`` appear on the queue.
#. ``status.txt`` in the recipe directory in ``$MAST_SCRATCH`` should show that ``perfect_opt1`` is queued.
#.  If you forgot some step above (like you forgot to create the submitlist file) and are running into strange problems, delete the PhononNebTest... folder from ``$MAST_SCRATCH`` and start again from the beginning of this section.
#.  The ``$MAST_CONTROL`` folder gives you error messages and other information. See :doc:`Running MAST <5_0_runningmast>` for tips.



*************************
Unit testing
*************************

To run unit tests and verify that the MAST code is sound, go to the test directory in your MAST installation path (e.g. <python installation path>/lib/python2.7/site-packages/MAST/test) and run the command ::

    nosetests -v --nocapture

The ``nocapture`` option allows print statements.
The ``verbose`` option gives verbose results.

The development team may have designated some tests to be skipped. However, any errors should be reported to the development team.
