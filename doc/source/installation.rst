============
Installation
============
#. Copy the top MAST directory into your user directory. We will refer to the top MAST directory as ``//home/user/topmast`` in this documentation.

#. Set the MAST environment variables. The script ``initialize.py`` in ``//home/user/topmast`` will give a set of default values to paste into your terminal window or setup profile, such as ``//home/user/.bashrc``. If you put the commands in your profile, you will have to either run ``source ~/.bashrc`` or log out and log back in. VASP_PSP_DIR is not defaulted.

    #. MAST_INSTALL_PATH: This variable should be set to the installation directory, for example::
    
        export MAST_INSTALL_PATH=//home/user/topmast

    #. MAST_RECIPE_PATH: Initially, this variable should be set to the existing recipe_templates folder. As you develop additional recipes, you may want to change this variable to reflect your customized organization of recipes.::
        
        export MAST_RECIPE_PATH=//home/user/topmast/recipe_templates

    #. MAST_SCRATCH: This variable may be set to any directory. Default behavior is to generate ingredients under this MAST_SCRATCH directory, unless the match_scratch keyword is specified with an overriding path in the input file.

    #. PYTHONPATH: If this environment variable already exists, the installation directory should be appended. Otherwise, this variable can be set to the installation directory. Assuming PYTHONPATH already has some value  (use ``env`` to see a list of environment variables):: 
        
        export PYTHONPATH=$PYTHONPATH://home/user/topmast
        
    #. VASP_PSP_DIR: This variable is necessary if VASP and VASP pseudopotential files are being used. See the documentation for the Materials Project's **pymatgen** code. The VASP_PSP_DIR should be set to a path which contains folder such as POT_GGA_PAW_PBE (for functional PBE, or mast_xc PBE in :doc:`Ingredients <ingredients>`) or POT_GGA_PAW_PW91 (for functional PW91).

    #. PATH: This variable should be appended with the MAST bin directory, for example::
    
        export PATH=$PATH://home/user/topmast/bin




#. Copy the appropriate example queue and script files for your platform from section :ref:`platforms`. **MAST team, we need a special test for these so that someone can run them and see if they work.**

.. _platforms:

----------------
Platform Support
----------------
Queue and submission script commands are in ``//home/user/topmast/submit`` and may need to be heavily modified depending on the platform used. 
To customize the queue submission behavior, copy the two _example.py files into new .py files, removing "_example" and overwriting the default files.
 
The out-of-the-box PBS submission script is built using the following input file keywords (see :doc:`Ingredients <ingredients>`):

* mast_processors or a combination of mast_ppn and mast_nodes
* mast_queue
* mast_exec
* mast_walltime
* mast_memory
* the ingredient name
