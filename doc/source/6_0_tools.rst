####################################
MAST post-processing utilities
####################################

******************************
Defect formation energy
******************************

The defect formation energy tool goes through the output of finished recipes in $MAST_ARCHIVE and calculates defect formation energies. It is found in $MAST_INSTALL_PATH/tools. 

The defect formation energy tool will create a ``<recipe directory>_dfe_results`` directory in the directory from which it is called.

To run without prompts::

    python $MAST_INSTALL_PATH/tools/defect_formation_energy <DFT bandgap> <experimental bandgap>

where DFT bandgap is a float for an LDA or GGA bandgap, and experimental bandgap is a float for an experimental or more accurate hybrid calculation bandgap.

To run with prompts::

    python $MAST_INSTALL_PATH/tools/defect_formation_energy prompt

*  Select the desired recipe
*  Follow the prompts for chemical potential conditions, band gap energy levels, and band gaps for adjustment

The two-column printout is Fermi energy on the left, and defect formation energy on the right.

*************************
Diffusion coefficient
*************************
**Update this section with new diffusion coefficient help text**

**************************
Defect finder
**************************

The defect finder takes a POSCAR file and finds vacancies and interstitials.
The defect finder currently exists in a separate repository.
You may test it online at materialshub.org > Resources > Tools > Defect Finder


