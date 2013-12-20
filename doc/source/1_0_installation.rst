============
Installation
============

------------------------
Installation on bardeen
------------------------
#. MAST is installed in ``//share/apps/MAST``.
#. Set the MAST environment variables. Add the following lines to your setup profile, such as ``//home/username/.bashrc``, where ``username`` is your username. Replace all instances of ``//home/username`` with your actual username, like ``//home/janedoe``. The environment variables are::
    
    export MAST_INSTALL_PATH=//share/apps/MAST
    export MAST_RECIPE_PATH=//home/username/MAST/recipe_templates
    export MAST_SCRATCH=//home/username/MAST/SCRATCH
    export MAST_ARCHIVE=//home/username/MAST/ARCHIVE
    export MAST_CONTROL=//home/username/MAST/CONTROL
    export PYTHONPATH=$PYTHONPATH://share/apps/MAST
    export VASP_PSP_DIR=//share/apps/MAST/vasp_pps
    export PATH=$PATH://share/apps/MAST/bin
    export PATH=//share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin:$PATH

#.   Explanation of environment variables:

    *  MAST_INSTALL_PATH: This variable should be set to the installation directory. ::

        export MAST_INSTALL_PATH=//share/apps/MAST

    *  MAST_RECIPE_PATH: MAST looks for recipe templates in this folder. You may want to copy recipes from the ``$MAST_INSTALL_PATH/recipe_templates`` directory into this folder and modify them. ::
    
        export MAST_RECIPE_PATH=//home/username/MAST/recipe_templates

    *  MAST_SCRATCH: This variable may be set to any directory. MAST will look for recipes in this directory. ::
    
        export MAST_SCRATCH=//home/username/MAST/SCRATCH

    *  MAST_ARCHIVE: This variable may be set to any directory. MAST will move completed recipes from ``$MAST_SCRATCH`` into this directory. ::
    
        export MAST_ARCHIVE=//home/username/MAST/ARCHIVE

    *  MAST_CONTROL: This variable may be set to any directory. MAST monitor log files, MAST monitor error files, and other MAST monitor output will be written to this directory. ::
    
        export MAST_CONTROL=//home/username/MAST/CONTROL

    *  PYTHONPATH: If this environment variable already exists, the MAST installation directory should be appended. Otherwise, this variable can be set to the installation directory. Assuming PYTHONPATH already has some value (use env to see a list of environment variables)::
    
        export PYTHONPATH=$PYTHONPATH://share/apps/MAST

    *  VASP_PSP_DIR: This variable is necessary if VASP and VASP pseudopotential files are being used. See the documentation for the `Materials Project's <http://materialsproject.org>` `pymatgen <http://pymatgen.org>` code. The VASP_PSP_DIR should be set to a path which contains folder such as POT_GGA_PAW_PBE (for functional PBE, or mast_xc PBE in Ingredients) or POT_GGA_PAW_PW91 (for functional PW91). ::
    
        export VASP_PSP_DIR=//share/apps/MAST/vasp_pps

    *  PATH: This variable should be appended with the ``$MAST_INSTALL_PATH/bin`` directory, for example::
    
        export PATH=$PATH://share/apps/MAST/bin:PATH

    
        *  Also, make sure that the correct version of python is defaulted to be used first. If you already use python for something else and this next line interferes with your other python calls (for example, you routinely use Python 2.4.3 instead and your other programs break if called from python 2.7.3), please see Tam. ::
        
        export PATH=//share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin:$PATH
    

        *  This python has pymatgen, numpy, and scipy in the appropriate libraries, which we need.
        *  Type ``which python`` and you should get: ``/share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python``
        *  Type ``which mast`` and you should get: ``/share/apps/MAST/bin/mast``

#.  Create all directories which do not yet exist (e.g., ``mkdir //home/username/MAST``, ``mkdir //home/username/MAST/recipe_templates``, ARCHIVE, CONTROL, and SCRATCH)
#.  Make an empty file at ``//home/username/MAST/CONTROL/submitlist``
#.  Log out of all bardeen terminals and log back in. (You may also run ``source ~/.bashrc``, but sometimes this doesn't quite set everything.)
#.  (There are some additional Platform Support steps which have already been taken: Queue and submission script commands are in $MAST_INSTALL_PATH/submit and may need to be heavily modified depending on the platform used. To customize the queue submission behavior, copy the appropriate files out of $MAST_INSTALL_PATH/submit/platforms and into $MAST_INSTALL_PATH/submit, omitting the platform name, and modify the new queue_commands.py, script_commands.py, and submit.sh accordingly. This has already been done on bardeen. No step here.)

---------------------------------
Test that MAST can run on bardeen
---------------------------------
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

#. Run mast once: ``nice -n 19 mast``
#. You should see a `mastmon` job appear on morganshort.
#. MAST should have detected that the first ingredient was ready to run, so when that process disappears, run mast again: ``nice -n 19 mast``
#. Now you should see ``perfect_opt1`` appear on the queue.
#. ``status.txt`` in the recipe directory in ``$MAST_SCRATCH`` should show that ``perfect_opt1`` is queued.
#. If you forgot some step above (like you forgot to create the submitlist file) and are running into strange problems, delete the PhononNebTest... folder from ``$MAST_SCRATCH`` and start again from the beginning of this section.
#. The ``$MAST_CONTROL`` folder gives you error messages and other information. See :doc:`Troubleshooting <5_0_troubleshooting>` for tips.

---------------------------------
Installation on another cluster
---------------------------------
1.  (On ACI/HPC, make sure you are using the compile node for all installation tasks. Use the submit node only to submit jobs.)
2.  Have the owner of //tmp/pip-build remove the directory if it exists (https://github.com/pypa/pip/issues/729
a.  cd //tmp
b.  rm -r pip-build
3.  Locate your version of python 2.7.3
a.  On platforms with .modules. it is probably something like .module load python. but get the correct version (.module avail. to see available modules). Type .which python. to make sure you have the right version, or .python --version.
i.  DLX has python 2.6.6 normally. .module load Python,. even though it is 2.7.3, has some difficulties installing pymatgen, possibly because of the way the module system works. Follow the .install python. directions instead.
b.  On bardeen it is //share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64
4.  If you do not have python, install it. The EPD/Canopy version is preferred because it includes numpy and scipy already
a.  https://www.enthought.com/downloads/
i.  version 2.7.5 is okay
b.  srun -u bash -i (on DLX, for interactive setup)
c.  bash ./canopy-1.0.3-rh5-64.sh
i.  Follow the prompts
d.  Add lines to your profile to make this your default python
i.  vi ~/.bashrc
ii. #EPD (Canopy) python
iii.    export PATH=//home/tma249/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin:$PATH
iv. Do not just use the .Canopy/bin. directory - python modules will not load properly
v.  Log out and log in
e.  Check your version of python: python --version
i.  This must be the correct version. If not, for all commands below which use .python,. give the full path to your version of python, e.g. //share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python
f.  Get setuptools (easy_install)
i.  wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
ii. python ez_setup.py
g.  Get pip
i.  curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
ii. python get-pip.py
h.  easy_install and pip are now wherever your installed python is.
i.  Check if easy_install and pip are available:
i.  which pip
ii. which easy_install
iii.    Example:
1.  [username@aci-service-2 ~]$ which pip
2.  //home/username/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin/pip
3.  [username@aci-service-2 ~]$ which easy_install
4.  //home/username/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin/easy_install
iv. pip must be version 1.3 or later (pip --version)
j.  If pip is not available and you are using the default version of python (not a local installation)
i.  You may need setuptools first:
1.  https://pypi.python.org/pypi/setuptools/0.9.8#installation-instructions
2.  wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
3.  python ez_setup.py --user
a.  Remember to use the correct version of python. Your actual python may be //home/<username>/bin/python-x.x.x/bin/python
b.  Or //share/apps/EPD...
ii. https://pypi.python.org/pypi/pip
iii.    http://www.pip-installer.org/en/latest/installing.html
iv. Option 1:
1.  curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
2.  python get-pip.py
v.  Option 2:
1.  curl -O https://pypi.python.org/packages/source/p/pip/pip-1.4.1.tar.gz
2.  nice -n 19 tar -xzvf pip-1.4.1.tar.gz
3.  cd pip.1.4.1
4.  python setup.py install --user
a.  Remember to use the correct version of python. Your actual python may be //home/<username>/bin/python-x.x.x/bin/python
b.  Or //share/apps/EPD...
vi. Now use the versions of easy_install and/or pip which are located in //home/<username>/.local/bin/
5.  Check if numpy is available:
a.  python (Use the correct version of python)
b.  import numpy
c.  If numpy is not available:
i.  Try pip installation. Depending on where pip is located, use the command:
1.  pip install --user numpy (Use the correct version of pip)
ii. Or use the command:
1.  //home/<username>/.local/bin/pip install --user numpy
iii.    If pip does not work, follow Quick install of numpy here:
iv. This will install Numpy without external library support.  It is a quick and easy way to install Numpy, and will suite you for the purposes of running MAST.
v.  Grab the most recent stable release of numpy
1.  http://www.scipy.org/install.html
vi. Untar with tar -zxvf numpy-<version>.tar.gz
vii.    cd numpy-<version>
viii.   Put the following in your command line:
1.  BLAS=None LAPACK=None ATLAS=None python setup.py config build install --prefix=<location where you want numpy installed, recommend $HOME/lib>
ix. Get something to drink, this.ll take about 5-10 minutes.
x.  Add to your .bashrc:
1.  NUMPY=<location you specified above>
2.  export PYTHONPATH=$NUMPY:$PYTHONPATH
xi. source $HOME/.bashrc
6.  Install pymatgen and custodian
a.  tma249@dlxlogin2-2 mast_installation]$ which pip
b.  //home/tma249/Canopy/appdata/canopy-1.0.3.1262.rh5-x86_64/bin/pip
c.  If .which easy_install. and .which pip. return the correct values, run the following commands.
i.  Otherwise, make sure you explicitly use the correct pip and easy_install, e.g. //home/username/.local/bin/pip and //home/username/.local/bin/easy_install or other such paths.
d.  Use the .--user. tag if you are not using the easy_install and pip from your own installation of python. Otherwise, you can omit this tag.
i.  nice -n 19 easy_install --user --upgrade distribute
1.  You MUST upgrade distribute, even if it is freshly installed. Just installing it will not work (8/9/13)
ii. nice -n 19 pip install --user pymatgen
iii.    nice -n 19 pip install --user custodian
e.  If pip does not work, try making your own temp directory.
i.  mkdir //home/<username>/tmp
ii. export TMPDIR=.//home/<username>/tmp.
iii.    Try running the pip commands again.
f.  If pymatgen fails to install, re-run steps (3.b.i) and (3.b.ii) again. Make sure that distribute has been upgraded.
7.  Remove any pip directory if it exists.
a.  cd //tmp
b.  rm -r pip-build
8.  Set up pymatgen VASP_PSP_DIR
a.  Locate the VASP pseudopotentials
i.  On bardeen, this is //share/apps/vasp_pseudopotentials
ii. On DLX it is //home/adozier/VASP
1.  On DLX, SKIP TO STEP 7.e
b.  Run pymatgen.s python setup tool
i.  This should now be wherever pymatgen was installed, either ~/.local/bin/potcar_setup.py if you installed it with --user, or wherever python is, otherwise.
ii. python .local/bin/potcar_setup.py or python potcar_setup.py or simply potcar_setup.py
iii.    (Remember to use the correct version of python, determined in step 2, e.g. //share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python .local/bin/potcar_setup.py)
c.  Example:
i.  Please enter full path where the POT_GGA_PAW_PBE, etc. subdirs are present. If you obtained the PSPs directly from VASP, this should typically be the directory that you untar the files to : //share/apps/vasp_pseudopotentials/paw
ii. Take the paw directory if you are using PAW. Do not take the top directory, or the GGA/LDA/etc folders will overwrite.
iii.    Please enter the fullpath of the where you want to create your pymatgen resources directory:
iv. //home/<username>/.local/vasp_pps
d.  Rename the folders under //home/<username>/.local/vasp_pps:
1.  rename the PBE folder POT_GGA_PAW_PBE to correspond to mast_xc pbe
2.  rename the GGA folder POT_GGA_PAW_PW91 to correspond to mast_xc pw91
a.  Add a line to your .bashrc file exporting the environment variable VASP_PSP_DIR to this VASP directory.
i.  On bardeen, it should look something like:
1.  export VASP_PSP_DIR=//home/<username>/.local/vasp_pps
ii. On DLX, use the directories already created:
1.  export VASP_PSP_DIR=//home/adozier/VASP/resources
iii.    or export VASP_PSP_DIR=<whichever path you used in the potcar_setup.py script>
iv. Remember to save your .bashrc file.
b.  Test the change:
i.  source ~/.bashrc
ii. cd $VASP_PSP_DIR
iii.    Make sure you are getting to the right directory, which has POT_GGA_POW_PBE etc. folders inside it.
9.  Make bin executables runnable:
a.  chmod -R a+x $MAST_INSTALL_PATH/bin
10. Modify submission details for your platform
a.  Go to $MAST_INSTALL_PATH/submit
b.  cp platforms/script_commands_<yourplatform>.py script_commands.py
c.  cp platforms/queue_commands_<yourplatform>.py queue_commands.py
d.  cp platforms/submit_<yourplatform>.sh submit.sh
11. Modify submit.sh as necessary for your platform.
a.  The submit.sh script should be set up to run mastmon.py on the shortest wallclock, fastest-turnaround queue on your system (e.g. a serial queue, morganshort, etc.)
b.  Examples of special modifications for submit.sh:
i.  ACI/HPC, add line: #SBATCH --partition=univ
ii. Bardeen, add a line to tell control where to run the monitor: #PBS -q morganshort
12. Modify script_commands.py as necessary for your platform.
a.  ACI/HPC: in script_commands.py, near line 95, add line: myscript.data.append("#SBATCH --partition=univ " + "\n")
b.  Bardeen: in script_commands.py near line 95 add line: myscript.data.append("#PBS -q " + mast_queue + "\n")
13. Modify queue_commands.py as necessary for your platform. 
14. Figure out the correct mast_exec calls for your system, to be used in input.inp. Examples are below.
a.  Bardeen: mast_exec //opt/mpiexec/bin/mpiexec //share/apps/bin/vasp5.2_par_opt1  (or any of the other vasp executables) 
b.  ACI/HPC: mast_exec //home/tma249/bin/vaspmpirun
i.  where vaspmpirun is this script (I put it in dlx.s //tmp/to_Henry):
ii. [tma249@dlxlogin2-2 bin]$ cat vaspmpirun
iii.    #!/bin/bash
iv. export PERL5LIB=/opt/moab/lib/perl5
v.  export MIC_LD_LIBRARY_PATH=/share/cluster/RHEL6.2/x86_64/apps/intel/ict/composer_xe_2013.0.079/compiler/lib/mic
vi. export LD_LIBRARY_PATH=/share/cluster/RHEL6.2/x86_64/apps/openmpi/1.6.2/lib:/share/cluster/RHEL6.2/x86_64/apps/intel/ict/composer_xe_2013.0.079/compiler/lib/intel64:/share/cluster/RHEL6.2/x86_64/apps/intel/ict/composer_xe_2013.0.079/mkl/lib/intel64
vii.    export INTEL_MKL_LIBS=/share/cluster/RHEL6.2/x86_64/apps/intel/ict/composer_xe_2013.0.079/mkl/lib/intel64
viii.   export QTLIB=/usr/lib64/qt-3.3/lib
ix. PATH=$PATH://home/tma249/bin://home/tma249/bin/convaspTest
x.  export PATH
xi. VaspPath=//home/adozier/VASP/vasp.5.2
xii.    export OMP_NUM_THREADS=1
xiii.   ulimit -s unlimited
xiv.    ulimit -l unlimited
xv. #mpirun $VaspPath/vasp
xvi.    //share/cluster/RHEL6.2/x86_64/apps/openmpi/1.6.2/bin/mpirun $VaspPath/vasp
15. Modify ~/.bashrc if necessary
a.  ACI/HPC, add line: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH://opt/intel/lib/intel64
16. To ensure recipes are created correctly, add python whitespace tab stops to your ~/.vimrc file:
a.  " VIM settings for python in a group below:
b.  set tabstop=4
c.  set shiftwidth=4
d.  set smarttab
e.  set expandtab
f.  set softtabstop=4
g.  set autoindent
17. Follow the environment variable setup in a similar fashion to Installation on bardeen
18. Follow the testing instructions from Test that MAST can run on bardeen




#. Copy the appropriate example queue and script files for your platform from section :ref:`platforms`. **MAST team, we need a special test for these so that someone can run them and see if they work.**

.. _platforms:

----------------
Platform Support
----------------
Queue and submission script commands are in ``//home/user/topmast/submit`` and may need to be heavily modified depending on the platform used. 
To customize the queue submission behavior, copy a queue_commands.py, script_commands.py and submit.sh from ``$MAST_INSTALL_PATH/submit/platforms`` to ``$MAST_INSTALL_PATH/submit/``. Remove the platform name from the file names.
 
The out-of-the-box PBS submission script is built using the following input file keywords (see :doc:`Ingredients <ingredients>`):

* mast_processors or a combination of mast_ppn and mast_nodes
* mast_queue
* mast_exec
* mast_walltime
* mast_memory
* the ingredient name
