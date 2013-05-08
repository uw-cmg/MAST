============
Tutorial
============
#. Follow the :doc:`installation instructions <installation>`. You will also need to set the VASP_PSP_DIR environment variable to a path containing a directory POT_GGA_PAW_PBE with the unzipped or tar.gz set of pseudopotentials (e.g. POTCAR.Al.gz). This setup is required by pymatgen, which MAST uses.
#. Create a new folder ``//home/user/mast_test``.
#. Copy the input file below as ``//home/user/mast_test/test.inp``. For information on the sections, see :doc:`Input File <inputfile>`. Change the two keywords below in ingredients_global. For more information, see :doc:`Ingredients <ingredients>`::

    # Test file
    $mast
    program VASP
    system_name MgAl
    $end

    $structure
    coord_type fractional
    begin coordinates
    Mg 0.000000 0.000000 0.000000
    Al 0.500000 0.500000 0.500000
    end

    begin lattice
    3.0 0.0 0.0
    0.0 3.0 0.0
    0.0 0.0 3.0
    end
    $end

    $ingredients
    begin ingredients_global
    mast_ppn <CHANGE TO THE NUMBER OF PROCESSORS PER NODE>
    mast_exec <CHANGE TO THE COMMAND TO RUN VASP>
    mast_kpoints 3x3x3 M
    mast_xc pbe
    end

    begin optimize
    encut 300
    ibrion 2
    isif 3
    end
    $end

    $recipe
    recipe_file recipe_test.txt
    $end

#. In your new folder, run ``mast -i test.inp`` to create your test job files and add the calculation directories to your scratch directory.
#. Look in your scratch directory (``echo $MAST_SCRATCH`` to see where this is) and make sure that empty test folders were created.
#. Run ``mast`` to start the scheduler and run your jobs. The first job (opt1) should complete, then forward its structure to the second job (opt2) to run.
#. When both jobs have completed, the main test folder should be moved into the ``archive`` folder under your scratch directory.

.. _nextsteps:

------------------
Next Steps
------------------
#. Sample recipe templates are provided in ``//home/user/topmast/recipe_templates``, but you will probably need to create your own. List all the ingredients at the top, like in a cooking recipe, and then list the Parent/Child logic.  See :doc:`Recipe <recipe>` for more information.

#. To run VASP several times with different parameters, use different ingredients in your recipe. Most VASP ingredients can be created by inheriting from the Optimize ingredient, just like the OptimizeLowMesh ingredient was created. Copy the optimizelowmesh.py file into a new file in the ``//home/user/topmast/MAST/ingredients`` directory. Open it up and change OptimizeLowMesh into a new class name. **We need a user-editable spot to put ingredients in, and some kind of helper to help copy ingredients, like copy optimizelowmesh.py into mynewoptimize.py with a new classname.** Then, in the input file, specify your different VASP parameters in the ``$ingredients`` section, for each ingredient. 

#. If you want the scheduler to run automatically, set up your crontab to run the command ``mast`` at a specified interval (for VASP runs, this should probably be at least 1 or 2 hours between mast runs).
