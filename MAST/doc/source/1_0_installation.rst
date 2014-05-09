#############
Installation
#############
 
*********************************
Installation
*********************************

===========================
Pre-steps 
===========================
Skip this step if you are on bardeen.

*  If you are on ACI/HPC, make sure you are using the compile node for all installation tasks. (aci-service-2 as of Dec. 2013) Use the submit node only to submit jobs.

*  Have the owner of //tmp/pip-build remove the directory if it exists (see `pip issue 729 <https://github.com/pypa/pip/issues/729>`_)::

    cd //tmp
    rm -r pip-build

================================
Verify your Python version
================================
Check your version of python: ``python --version``

If your version of python is not 2.7.3, try to locate an existing version of python 2.7.3.
Then, make sure that this version of python is defaulted to be used first. You may need to add a line to your user profile. Your user profile may be located in ``//home/username/.bashrc`` or a similar file.

For bardeen, the line you need to add is::
    
    export PATH=//share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin:$PATH

Then, log out and log back in. 

For platforms with the "module" system like stampede or DLX, check which modules are available (``module avail``) and add a line something like::

    module load python
    module load Python

Then, log out and log back in.

Type ``which python`` to make sure you have the right version, or ``python --version``.

If you already use python for something else and shifting python versions will interfere with other programs, for example, you routinely use Python 2.4.3 instead and your other programs break if called from python 2.7.3, please contact the development team. ::
        
If you do not have or cannot find Python 2.7.3, then you must install it. 

---------------------
Install python
---------------------
The EPD/Canopy version is preferred because it includes numpy and scipy already. Download this version from `EPD Free Canopy <https://www.enthought.com/downloads/>`_

*  Version 2.7.5 is okay
*  On DLX, go into interactive setup with the command ``srun -u bash -i``
*  ``bash ./canopy-1.0.3-rh5-64.sh``
*  Follow the prompts (use spacebar to scroll through the license file)

Add lines to your profile to make this python installation your default python::

    vi ~/.bashrc
    #EPD (Canopy) python
    export PATH=//home/<username>/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin:$PATH

*  Do not just use the .Canopy/bin. directory - python modules will not load properly
*  Log out and log in

Check your version of python: ``python --version``

The version given must be the correct version. If not, for all subsequent commands that say *python*, give the full path to your version of python, e.g. ``//share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python``
        
==============================================
Verify setuptools (easy_install) and pip
==============================================
Check if easy_install and pip are available:

    *  ``which pip``
    *  ``which easy_install``

Example::
    
    [username@aci-service-2 ~]$ which pip
    //home/username/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin/pip
    [username@aci-service-2 ~]$ which easy_install
    //home/username/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin/easy_install
    
pip must be version 1.3 or later (``pip --version``)

If either easy_install or pip is missing, install them as follows.

Get setuptools (easy_install)

    *  `setuptools <https://pypi.python.org/pypi/setuptools>`_
    *  ``wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py``
    *  ``python ez_setup.py`` if you are using your own locally-installed python
    *  ``python ez_setup.py --user`` if you are using a root-installed python

Get pip

    *  `pip <https://pypi.python.org/pypi/pip>`_
    *  ``curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py``
    *  ``python get-pip.py`` if you are using your own locally-installed python
    *  ``python get-pip.py --user`` if you are using a root-installed python

easy_install and pip should now be located either wherever your installed python is, or in the ``$HOME/.local/bin`` directory
Check their locations and the pip version again.

=========================================
Verify or install numpy and scipy
=========================================
Check if numpy and scipy available::

    python
    import numpy
    import scipy

If numpy and scipy are not available, we recommend that you go back and install a local version of python which already includes numpy and scipy.

Scipy is optional at this stage (used in the MAST defect finder).

--------------------------------------------
Install numpy (not recommended)
--------------------------------------------  
If numpy is not available, try pip installation::

    pip install --user numpy

(If you are using a user-installed pip with a root-installed python, use the command ``$HOME/.local/bin/pip`` instead of ``pip``.)

If pip does not work, follow Quick install of numpy here. This will install Numpy without external library support. It is a quick and easy way to install Numpy, and will suite you for the purposes of running MAST.

    *  Grab the most recent stable release of numpy from `<http://www.scipy.org/install.html>`_
    *  Untar with command ``tar -zxvf numpy-<version>.tar.gz``
    *  ``cd numpy-<version>``
    *  Put the following in your command line, all as one line::

        BLAS=None LAPACK=None ATLAS=None 
        python setup.py config build install 
        --prefix=<location where you want numpy installed, recommend $HOME/lib>

    *  Get something to drink; this'll take about 5-10 minutes.
    *  Add to your .bashrc::
        
        NUMPY=<location you specified above>
        export PYTHONPATH=$NUMPY:$PYTHONPATH

    *  source $HOME/.bashrc

============================================
Verify or install pymatgen and custodian
============================================   
Check if pymatgen and custodian available::

    python
    import pymatgen
    import custodian

If pymatgen and custodian are not available, install them.

--------------------------------
Install pymatgen and custodian
--------------------------------

Make sure you explicitly use the correct pip and easy_install, e.g. //home/username/.local/bin/pip and //home/username/.local/bin/easy_install or other such paths, corresponding to the correct version of python.

Use the ``--user`` tag if you are not using the easy_install and pip from your own installation of python. Otherwise, you can omit this tag.

Upgrade the *distribute* package. You **MUST** upgrade this package, even if it is freshly installed. (8/9/13) ::
    
    nice -n 19 easy_install --user --upgrade distribute

pip install pymatgen and custodian::

    nice -n 19 pip install --user pymatgen
    nice -n 19 pip install --user custodian

If the pymatgen installation does not work, failing with PyCifRW, install PyCifRW manually first, using the paths that correspond to your system (python line is all one line)::

    cd $HOME/.local/lib/python2.7/site-packages/setuptools-2.1-py2.7.egg

    python ./easy_install.py --user https://bitbucket.org/
        jamesrhester/pycifrw/downloads/PyCifRW-3.5-py2.7-linux-i686.egg

If pip does not work, try making your own temp directory. ::
            
    mkdir //home/<username>/tmp
    export TMPDIR=.//home/<username>/tmp.

Then try running the pip commands again.
            
Remove any pip directory if it exists. ::
    
    cd //tmp
    rm -r pip-build

======================================
Set up the pymatgen VASP_PSP_DIR
======================================
On DLX and bardeen, skip to the NEXT NUMBERED STEP

Locate the VASP pseudopotentials

*  On bardeen, this is ``//share/apps/vasp_pseudopotentials``
*  On DLX it is ``//home/adozier/VASP``
    
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
Get MAST
===============================
* Get the latest MAST package from the Python package index::
    
    pip install --upgrade --no-deps MAST

(No dependencies because we are assuming pymatgen and custodian have been properly installed as above. It is recommended to install them separately.)
    
* The TESTING LINK is: pip install --upgrade --no-deps MAST_tam_test -i https://testpypi.python.org/pypi

======================================
Set up the environment variables
======================================
The pip installation should set up a ``MAST`` directory in ``//home/username/MAST`` with several subdirectories.

The pip installation should then warn you with an ATTENTION flag of environment variables that must be set. 

You may copy and paste the environment variables from the terminal into your user profile. In the examples below, ``username`` should have been changed to your username.::
    export MAST_RECIPE_PATH=//home/username/MAST/recipe_templates
    export MAST_SCRATCH=//home/username/MAST/SCRATCH
    export MAST_ARCHIVE=//home/username/MAST/ARCHIVE
    export MAST_CONTROL=//home/username/MAST/CONTRO"
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

MAST_RECIPE_PATH: MAST looks for recipe templates in this folder. You will have been supplied with a few example templates, also corresponding to example input files in a newly-created ``//home/username/MAST/examples`` directory. ::
    
    export MAST_RECIPE_PATH=//home/username/MAST/recipe_templates

MAST_SCRATCH: This variable may be set to any directory. MAST will look for recipes in this directory. ::
    
    export MAST_SCRATCH=//home/username/MAST/SCRATCH

MAST_ARCHIVE: This variable may be set to any directory. MAST will move completed recipes from ``$MAST_SCRATCH`` into this directory. ::
    
    export MAST_ARCHIVE=//home/username/MAST/ARCHIVE

MAST_CONTROL: This variable may be set to any directory. MAST monitor log files, MAST monitor error files, and other MAST monitor output will be written to this directory. ::
    
    export MAST_CONTROL=//home/username/MAST/CONTROL

MAST_CONTROL also has several subfolders. If you move your $MAST_CONTROL to a different path, please copy the subfolders with it.

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
#.  Make a test directory, like ``//home/username/MAST/test``
#.  Copy the test input file to your test folder (all one line)::

        cp //share/apps/MAST/test/phononreorgtest/
        phonon_with_neb.inp //home/username/MAST/test/test.inp

#.  Go to your test directory, ``cd //home/username/MAST/test``
#.  Modify the test.inp file with the correct ``mast_exec``, ``mast_ppn``, ``mast_queue``, and other settings described in :ref:`platforms`
#.  Try to parse the input file, entering the following command as one line::

        nice -n 19 mast -i test.inp 

    *  The .nice -n 19. keeps this command low priority, since it is being run on the headnode (but it is not too intensive).
    *  The -i signals to MAST that it is processing an input file.
#. Your ``//home/username/MAST/SCRATCH`` directory should now have a folder with a very long name in it (recipe directory), which contains several subfolders (ingredient directories).
#. Go to that long recipe directory. (PhononNebTest...)

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
#.  You should see a `mastmon` job appear on morganshort.
#.  MAST should have detected that the first ingredient was ready to run, so when that process disappears, run mast again: ``nice -n 19 mast``
#.  Now you should see ``perfect_opt1`` appear on the queue.
#. ``status.txt`` in the recipe directory in ``$MAST_SCRATCH`` should show that ``perfect_opt1`` is queued.
#.  If you forgot some step above (like you forgot to create the submitlist file) and are running into strange problems, delete the PhononNebTest... folder from ``$MAST_SCRATCH`` and start again from the beginning of this section.
#.  The ``$MAST_CONTROL`` folder gives you error messages and other information. See :doc:`Running MAST <5_0_runningmast>` for tips.



*************************
Unit testing
*************************

To run unit tests and verify that the MAST code is sound, go to
``$MAST_INSTALL_PATH/test`` and run the command ::

    nosetests -v --nocapture

Or, optionally, run the command ::

    nosetests -v --nocapture

The ``nocapture`` option allows print statements.
The ``verbose`` option gives verbose results.

The development team may have designated some tests to be skipped. However, any errors should be reported to the development team.
