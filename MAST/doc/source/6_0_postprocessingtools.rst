####################################
MAST post-processing utilities
####################################

These utilities are meant to be used as part of a MAST workflow.
See example files in ``$HOME/MAST/examples`` or wherever you may have moved the initally-created ``$HOME/MAST/examples`` folder for examples on how to use them.

These utilities should have been copied into your bin or .local/bin directory (see :doc:`1_0_installation`).

***********************************************************************
Defect formation energy with finite-size scaling
***********************************************************************

Initially determining the sizes for finite-size scaling is covered in :doc:`3_1_1_structure` with the utility ``mast_finite_size_scaling_sizes``.

The ``mast_madelung_utility`` utility runs as the last ingredient in a finite-size scaling defect workflow (see ``$HOME/MAST/examples/finite_size_scaling.inp``).

Run the utility as ``mast_madelung_utility``. All inputs are derived from the recipe-local ``input.inp`` file in the recipe directory.

*  The utility should generate a series of tables and .png plots that display the finite-size-scaling-corrected and original defect formation energies for different chemical potentials.

* :doc:`3_1_7_chemicalpotentials` of the input file should be set in order for the utility to work.

***********************************************************************
Defect formation energy versus Fermi energy
***********************************************************************

The ``mast_defect_formation_energy`` tool plots defect formation energy versus Fermi energy. 

The defect formation energy tool is intended to be run as another ingredient folder in the recipe directory.

If you do not have such an ingredient in the recipe directory, you may manually create the ingredient folder and give it a ``dfe_input.txt`` file.

The ``dfe_input.txt`` file for a manually-created or embedded workflow ingredient (see ``//home/<username>/MAST/example/defect_formation_energy.inp``) should contain the following information::

    dfe_label1=perfect_label defected_label
    dfe_label2=perfect_label defected_label
    dfe_label3=perfect_label defected_label
    (etc. for more defects)
    bandgap_lda_or_gga=<float>
    bandgap_hse_or_expt=<float>
    plot_threshold <float>: Plotting threshold value

*  <perfect_label> and <defected_label> are the ingredient names of the perfected and corresponding defected cells.

*  bandgap_lda_or_gga should be a float value indicating a DFT-calculated bandgap, usually expected to be underestimated.

*  bandgap_hse_or_expt should be a float value indicating an experimental or more accurate bandgap, e.g. from a hybrid calculation.

*  plot_threshold should be a float value indicating the threshold for transitions.

*  In addition, :doc:`3_1_7_chemicalpotentials` should exist in the ``input.inp`` input file inside the recipe directory.

Run the utility as::

    mast_defect_formation_energy dfe_input.txt

A directory named ``dfe_results`` should be created within the ingredient directory. Inside that directory:

*  The two-column printout for each chemical potential-labeled text file contains Fermi energy on the left, and defect formation energy on the right.

*  The ``dfe.txt`` printout contains defect formation energy information for each charge state.

*************************
Diffusion coefficient
*************************

The ``mast_diffusion_coefficient`` diffusion coefficient calculation tool supports the following models:

*  FCC five-frequency model equation from R. E. Howard and J. R. Manning, Physical Review, Vol. 154, 1967.

*  FCC concentrated fourteen-frequency model equation from Bocquet J.-L. and Le Claire A. D.

*  HCP eight-frequency model equation from P. B. Ghate, Physical Review, Vol. 133, 1963.

The tool is designed to be used as a separate ingredient within the recipe directory. See ``$HOME/MAST/examples/neb_with_phonons.inp`` for an example input fileof a full workflow.

If the ingredient was not created within the workflow, an ingredient directory may be manually created for the tool.

The tool will use an input text file like ``diffcoeff_input.txt``, which should contain the following lines. The order of the lines does not matter.

*  Names of the directories of energies and attempt rates, which are specified with respect to different frequencies for the model:
    
    *  **E** and **v** means energy and attempt rate, respectively. (There is no support for other characters such as w).

    *  For 5-freq, **E0 through E4** should be used to specify the relations with certain directories

    *  For 8-freq, **Ea, Eb, Ec, EX, Eap (p means prime), Ebp, Ecp, and EXp** should be used. Note they are all case sensitive and should be exactly the same as written here.

    *  Generally speaking, each keyword (Exx or vxx) is followed by two ingredient names. 
    
        *  The first name indicates the ingredient name corresponding to the configuration of the starting point of NEB.
        
        *  The second name indicates the ingredient name corresponding to the configuration of the saddle point of the NEB.
        
        *  This order should not be changed.

        *  For each name, the utility will expect two files to be present within the ingredient diretory of the diffusion coefficient tool:
        
            * <ingredient_name>_OUTCAR

            * <ingredient_name>_OSZICAR

            * If you are manually creating a diffusion coefficient tool ingredient, you will have to manually copy files from each of the completed ingredients specified.

    *  The user can also type only one single float behind the keyword, and the code will then not refer to the directory for the related energy or attempting rate, but simply use the data given.


*  **type** means which frequency model to choose. Either ``5`` or ``fcc`` tells the code that the five-frequency model should be applied, while either ``8`` or ``hcp`` tell the code that the eight-frequency model should be applied.


*  **HVf** means the formation energy of the vacancy

    * Either 1 float or two ingredient names are expected after this keyword.

    * If ingredient names are used, in the order <perfect_ingredient> <defected_ingredient>, then the utility will expect two energy files to be present in the utility's ingredient directory:

        * <perfect_ingredient>_OSZICAR

        * <defected_ingredient>_OSZICAR

        * Charged defects are not currently supported.

*  **HB** means the binding energy, and is only applicable for the 8-frequency model.

    * Either 1 float or four ingredient names are expected after this keyword.
    
    * If ingredient names are used:

        * Use the order <perfect ingredient> <vacancy and substitution> <substitution only> <vacancy only>

        * Supply an <ingredient_name>_OSZICAR file in the utility's ingredient directory.
    
*  **lattice** indicates the ingredient name for the ingredient in which to find a lattice file.

    *  This ingredient typically corresponds to an undefected supercell. 
    
    *  The utility expects to find a <lattice_ingredient_name>_POSCAR file inside the diffusion coefficient utility ingredient directory.

*  **plotdisplay** indicates whether to use matplotlib.pyplot in order to create a plot, or whether to skip plotting. 

    *  Use "plotdisplay none" to skip plotting

    *  Omit this keyword to use a default display

    *  Use "plotdiplay tkagg" etc. or another display string to specify a matplotlib display.

Run as ``mast_diffusion_coefficient -i <input>``

