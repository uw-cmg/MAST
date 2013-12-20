#############
Installation
#############

.. _installation-on-bardeen:

************************
Installation on bardeen
************************

========================
1. Environment variables
========================

MAST is installed in ``//share/apps/MAST``.

Set the MAST environment variables. Add the following lines to your setup profile, such as ``//home/username/.bashrc``, where ``username`` is your username. Replace all instances of ``//home/username`` with your actual username, like ``//home/janedoe``. 

The environment variables are::
   
    export MAST_INSTALL_PATH=//share/apps/MAST
    export MAST_RECIPE_PATH=//home/username/MAST/recipe_templates
    export MAST_SCRATCH=//home/username/MAST/SCRATCH
    export MAST_ARCHIVE=//home/username/MAST/ARCHIVE
    export MAST_CONTROL=//home/username/MAST/CONTROL
    export PYTHONPATH=$PYTHONPATH://share/apps/MAST
    export VASP_PSP_DIR=//share/apps/MAST/vasp_pps
    export PATH=$PATH://share/apps/MAST/bin
    export PATH=//share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin:$PATH

An explanation of each variable appears below.

MAST_INSTALL_PATH: This variable should be set to the installation directory. ::

    export MAST_INSTALL_PATH=//share/apps/MAST

MAST_RECIPE_PATH: MAST looks for recipe templates in this folder. You may want to copy recipes from the ``$MAST_INSTALL_PATH/recipe_templates`` directory into this folder and modify them. ::
    
    export MAST_RECIPE_PATH=//home/username/MAST/recipe_templates

MAST_SCRATCH: This variable may be set to any directory. MAST will look for recipes in this directory. ::
    
    export MAST_SCRATCH=//home/username/MAST/SCRATCH

MAST_ARCHIVE: This variable may be set to any directory. MAST will move completed recipes from ``$MAST_SCRATCH`` into this directory. ::
    
    export MAST_ARCHIVE=//home/username/MAST/ARCHIVE

MAST_CONTROL: This variable may be set to any directory. MAST monitor log files, MAST monitor error files, and other MAST monitor output will be written to this directory. ::
    
    export MAST_CONTROL=//home/username/MAST/CONTROL

PYTHONPATH: If this environment variable already exists, the MAST installation directory should be appended. Otherwise, this variable can be set to the installation directory. Assuming PYTHONPATH already has some value (use env to see a list of environment variables)::
    
    export PYTHONPATH=$PYTHONPATH://share/apps/MAST

VASP_PSP_DIR: This variable is necessary if VASP and VASP pseudopotential files are being used. See the documentation for the :ref:`Materials Project's <http://materialsproject.org>`_ :ref:`pymatgen <http://pymatgen.org>`_ code. The VASP_PSP_DIR should be set to a path which contains folder such as POT_GGA_PAW_PBE (for functional PBE, or mast_xc PBE in Ingredients) or POT_GGA_PAW_PW91 (for functional PW91). ::
    
    export VASP_PSP_DIR=//share/apps/MAST/vasp_pps

PATH: This variable should be appended with the ``$MAST_INSTALL_PATH/bin`` directory, for example::
    
    export PATH=$PATH://share/apps/MAST/bin:PATH

====================
2. Python version
====================

Make sure that the correct version of python is defaulted to be used first. If you already use python for something else and this next line interferes with your other python calls (for example, you routinely use Python 2.4.3 instead and your other programs break if called from python 2.7.3), please see Tam. ::
        
    export PATH=//share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin:$PATH
    
This version of python has pymatgen, numpy, and scipy in the appropriate versions.
    
*  Type ``which python`` and you should get: ``/share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python``
*  Type ``which mast`` and you should get: ``/share/apps/MAST/bin/mast``

=====================
3. Additional setup
=====================
#.  Create all directories which do not yet exist (e.g., ``mkdir //home/username/MAST``, ``mkdir //home/username/MAST/recipe_templates``, ARCHIVE, CONTROL, and SCRATCH)
#.  Make an empty file at ``//home/username/MAST/CONTROL/submitlist``
#.  Log out of all bardeen terminals and log back in. (You may also run ``source ~/.bashrc``, but sometimes this doesn't quite set everything.)
#.  (There are some additional Platform Support steps which have already been taken: Queue and submission script commands are in $MAST_INSTALL_PATH/submit and may need to be heavily modified depending on the platform used. To customize the queue submission behavior, copy the appropriate files out of $MAST_INSTALL_PATH/submit/platforms and into $MAST_INSTALL_PATH/submit, omitting the platform name, and modify the new queue_commands.py, script_commands.py, and submit.sh accordingly. This has already been done on bardeen. No step here.)

.. _test-on-bardeen:

*********************************
Test that MAST can run on bardeen
*********************************
#.  Copy the test recipe template to your recipe_templates folder::

        cp //share/apps/MAST/recipe_templates/phonon_test_neb.txt //home/username/MAST/recipe_templates/.

#.  Make a test directory, like ``//home/username/MAST/test``
#.  Copy the test input file to your test folder::

        cp //share/apps/MAST/test/phononreorgtest/phonon_with_neb.inp //home/username/MAST/test/test.inp

#.  Go to your test directory, ``cd //home/username/MAST/test``
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
#.  The ``$MAST_CONTROL`` folder gives you error messages and other information. See :doc:`Troubleshooting <5_0_troubleshooting>` for tips.
 
*********************************
Installation on another cluster
*********************************

===========================
0. Pre-steps
===========================
*  If you are on ACI/HPC, make sure you are using the compile node for all installation tasks. (aci-service-2 as of Dec. 2013) Use the submit node only to submit jobs.

*  Have the owner of //tmp/pip-build remove the directory if it exists (see :ref:`pip issue 729 <https://github.com/pypa/pip/issues/729>`_)::

    cd //tmp
    rm -r pip-build

================================
1. Verify your Python version
================================
Check your version of python: ``python --version``
If your version of python is not 2.7.3, try to locate an existing version of python 2.7.3.

*  On platforms with modules, it is probably something like ``module load python``, but get the correct version (``module avail`` to see available modules). Type ``which python`` to make sure you have the right version, or ``python --version``.

*  DLX has python 2.6.6 normally. ``module load Python``, even though it is 2.7.3, has some difficulties installing pymatgen, possibly because of the way the module system works. Follow the ``install python`` directions instead.

*  On bardeen it is //share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64

If you do not have python 2.7.3, install it. 

---------------------
Installing python
---------------------
The EPD/Canopy version is preferred because it includes numpy and scipy already. Download this version from here: :ref:`EPD Free Canopy <https://www.enthought.com/downloads/>`_

*  Version 2.7.5 is okay
*  On DLX, go into interactive setup with the command ``srun -u bash -i``
*  ``bash ./canopy-1.0.3-rh5-64.sh``
*  Follow the prompts (use spacebar to scroll through the license file)

Add lines to your profile to make this python installation your default python::

    vi ~/.bashrc
    #EPD (Canopy) python
    export PATH=//home/tma249/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin:$PATH

*  Do not just use the .Canopy/bin. directory - python modules will not load properly
*  Log out and log in

Check your version of python: ``python --version``

The version given must be the correct version. If not, for all subsequent commands that say *python*, give the full path to your version of python, e.g. ``//share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python``
        
==============================================
2. Verify setuptools (easy_install) and pip
==============================================
Check if easy_install and pip are available:

    *  ``which pip``
    *  ``which easy_install``

Example::
    
    [username@aci-service-2 ~]$ which pip
    //home/username/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin/pip
    [username@aci-service-2 ~]$ which easy_install
    //home/username/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin/easy_install
    
pip must be version 1.3 or later (pip --version)

If either easy_install or pip is missing, install them as follows.

Get setuptools (easy_install)

    *  :ref:`setuptools <https://pypi.python.org/pypi/setuptools>`_
    *  ``wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py``
    *  ``python ez_setup.py``

Get pip

    *  :ref:`pip <https://pypi.python.org/pypi/pip>`_
    *  ``curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py``
    *  ``python get-pip.py``

easy_install and pip should now be located wherever your installed python is.
Check their locations and the pip version again.

=================================
3. Verify numpy and scipy
=================================
Check if numpy and scipy available::

    python
    import numpy
    import scipy

  
If numpy is not available, try pip installation::

    pip install --user numpy

(Use the pip in the bin directory of the correct version of python)

If pip does not work, follow Quick install of numpy here. This will install Numpy without external library support. It is a quick and easy way to install Numpy, and will suite you for the purposes of running MAST.

    *  Grab the most recent stable release of numpy from :ref:`<http://www.scipy.org/install.html>`
    *  Untar with command ``tar -zxvf numpy-<version>.tar.gz``
    *  ``cd numpy-<version>``
    *  Put the following in your command line::

        BLAS=None LAPACK=None ATLAS=None python setup.py config build install --prefix=<location where you want numpy installed, recommend $HOME/lib>

    *  Get something to drink; this'll take about 5-10 minutes.
    *  Add to your .bashrc::
        
        NUMPY=<location you specified above>
        export PYTHONPATH=$NUMPY:$PYTHONPATH

    *  source $HOME/.bashrc

==================================
4. Install pymatgen and custodian
==================================    
Make sure you explicitly use the correct pip and easy_install, e.g. //home/username/.local/bin/pip and //home/username/.local/bin/easy_install or other such paths, corresponding to the correct version of python

Use the ``--user`` tag if you are not using the easy_install and pip from your own installation of python. Otherwise, you can omit this tag.

Upgrade the *distribute* package. You **MUST** upgrade this package, even if it is freshly installed. (8/9/13) ::
    
    nice -n 19 easy_install --user --upgrade distribute

pip install pymatgen and custodian::

    nice -n 19 pip install --user pymatgen
    nice -n 19 pip install --user custodian

If pip does not work, try making your own temp directory. ::
            
    mkdir //home/<username>/tmp
    export TMPDIR=.//home/<username>/tmp.

Then try running the pip commands again.
            
Remove any pip directory if it exists. ::
    
    cd //tmp
    rm -r pip-build

======================================
5. Set up the pymatgen VASP_PSP_DIR
======================================

Locate the VASP pseudopotentials

*  On bardeen, this is ``//share/apps/vasp_pseudopotentials``
*  On DLX it is ``//home/adozier/VASP``

    *  On DLX, SKIP TO THE NEXT NUMBERED STEP
    
Run pymatgen's python setup tool. This tool should be located wherever pymatgen was installed, either ``~/.local/bin/potcar_setup.py`` if you installed it with ``--user``, or wherever python is, otherwise. ::

    python .local/bin/potcar_setup.py or python potcar_setup.py or simply potcar_setup.py
        
(Remember to use the correct version of python, determined in step 2, e.g. //share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python .local/bin/potcar_setup.py)

Take the paw directory if you are using PAW. Do not take the top directory, or the GGA/LDA/etc folders will overwrite.

Example of running the python setup tool::
        
    Please enter full path where the POT_GGA_PAW_PBE, etc. subdirs are present. If you obtained the PSPs directly from VASP, this should typically be the directory that you untar the files to : //share/apps/vasp_pseudopotentials/paw
    Please enter the fullpath of the where you want to create your pymatgen resources directory:
    //home/<username>/.local/vasp_pps

Rename the folders under ``//home/<username>/.local/vasp_pps``:
    
*  Rename the PBE folder POT_GGA_PAW_PBE to correspond to mast_xc pbe
*  Rename the GGA folder POT_GGA_PAW_PW91 to correspond to mast_xc pw91

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
6. Get MAST
===============================

Get the MAST.tar.gz file from MaterialsHub.org and untar it::

    nice -n 19 tar -xzvf MAST.tar.gz


(or run this command over interactive submission, which is better)

Make the bin executables runnable. Supposing the uncompressed path was //home/username/MAST, then run the following command::

    chmod -R a+x //home/username/MAST/bin

Modify the submission details for your platform::

    cd //home/username/MAST/submit
    cp platforms/script_commands_<yourplatform>.py script_commands.py
    cp platforms/queue_commands_<yourplatform>.py queue_commands.py
    cp platforms/submit_<yourplatform>.sh submit.sh

Modify submit.sh as necessary for your platform.

*  The submit.sh script should be set up to run mastmon.py on the shortest wallclock, fastest-turnaround queue on your system (e.g. a serial queue, morganshort, etc.)
*  Examples of special modifications for submit.sh:
        
    *  ACI/HPC, add line: ``#SBATCH --partition=univ``
    *  Bardeen, add a line to tell control where to run the monitor: ``#PBS -q morganshort``

Modify script_commands.py as necessary for your platform.

*  ACI/HPC: in script_commands.py, near line 95, add line: ``myscript.data.append("#SBATCH --partition=univ " + "\n")``
*  Bardeen: in script_commands.py near line 95 add line: ``myscript.data.append("#PBS -q " + mast_queue + "\n")``

Modify queue_commands.py as necessary for your platform. (On DLX, ACI, and bardeen, no modification should be necessary.)

Figure out the correct mast_exec calls for your system, to be used in the :doc:`Input File<3_0_inputfile>`. Examples are below.

*  Bardeen: ``mast_exec //opt/mpiexec/bin/mpiexec //share/apps/bin/vasp5.2_par_opt1``  (or any of the other vasp executables) 
*  DLX: ``mast_exec //home/username/bin/vaspmpirun``, where vaspmpirun is the following script::

    #!/bin/bash
    export PERL5LIB=/opt/moab/lib/perl5
    export MIC_LD_LIBRARY_PATH=/share/cluster/RHEL6.2/x86_64/apps/intel/ict/composer_xe_2013.0.079/compiler/lib/mic
    export LD_LIBRARY_PATH=/share/cluster/RHEL6.2/x86_64/apps/openmpi/1.6.2/lib:/share/cluster/RHEL6.2/x86_64/apps/intel/ict/composer_xe_2013.0.079/compiler/lib/intel64:/share/cluster/RHEL6.2/x86_64/apps/intel/ict/composer_xe_2013.0.079/mkl/lib/intel64
    export INTEL_MKL_LIBS=/share/cluster/RHEL6.2/x86_64/apps/intel/ict/composer_xe_2013.0.079/mkl/lib/intel64
    export QTLIB=/usr/lib64/qt-3.3/lib
    PATH=$PATH:$HOME/bin:$HOME/bin/convaspTest
    export PATH
    VaspPath=//home/adozier/VASP/vasp.5.2
    export OMP_NUM_THREADS=1
    ulimit -s unlimited
    ulimit -l unlimited
    #mpirun $VaspPath/vasp
    //share/cluster/RHEL6.2/x86_64/apps/openmpi/1.6.2/bin/mpirun $VaspPath/vasp

Modify ~/.bashrc if necessary
    
*  ACI/HPC, add line: ``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH://opt/intel/lib/intel64``

To ensure recipes are created correctly, add python whitespace tab stops to your ~/.vimrc file::
    
    " VIM settings for python in a group below:
    set tabstop=4
    set shiftwidth=4
    set smarttab
    set expandtab
    set softtabstop=4
    set autoindent

Follow the environment variable setup in a similar fashion to :ref:`installation-on-bardeen`
Follow the testing instructions from :ref:`test-on-bardeen`



.. _platforms:

****************
Platform Support
****************
Queue and submission script commands are in ``//home/user/topmast/submit`` and may need to be heavily modified depending on the platform used. 
To customize the queue submission behavior, copy a queue_commands.py, script_commands.py and submit.sh from ``$MAST_INSTALL_PATH/submit/platforms`` to ``$MAST_INSTALL_PATH/submit/``. Remove the platform name from the file names.
 
The out-of-the-box PBS submission script is built using the following input file keywords (see :doc:`Ingredients <2_0_ingredients>`):

* mast_processors or a combination of mast_ppn and mast_nodes
* mast_queue
* mast_exec
* mast_walltime
* mast_memory
* the ingredient name
