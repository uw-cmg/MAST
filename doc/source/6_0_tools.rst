.. _6_0_tools:

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
Usage of diffusion coefficient calculation tool code
1.  This code currently supports 5(fcc) and 8(hcp) frequency models.
2.  The code currently will work in the same directory with other MAST generated folders (``neb_vac*``, ``phonon_vac*``, etc.)
3.  Type ``$MAST_INSTALL_PATH/MAST/utility/diff.py -i <input>`` to run.
4.  The input file should contain the following lines, naming the directories of energies and attempt rates which are specified with respect to different frequencies for the model.
a)  The order of different lines does not matter.
b)  There can be as many .\n. between lines or as many spaces between words, and they will not affect the code. 
c)  The keyword at the beginning of each line matters:
i.  .type. means which frequency model to choose. Either .5. or .fcc. tells the code 5-freq is applied, while either .8. or .hcp. the 8-freq.
ii. .E*. and .v*. means energy and attempting rate, respectively. (Currently doesn.t support other characters such as .w.).
iii.    For 5-freq, .E0.~.E4. should be used to specify their relations with certain directories; for 8-freq, .Ea.,.Eb.,.Ec.,.EX.,.Eap.(.p. means .prime.),.Ebp.,.Ecp. and .EXp. should be used. Note they are all case sensitive and should be exactly the same as 
iv. Generally speaking, each keyword (.E*. or .v*.) is followed by two words. The first indicates the configuration of the starting/end point of NEB and the second represents the saddle point. This order should not be changed.
v.  The user can also type only one single float behind the keyword, and the code will then not refer to the directory for the related energy or attempting rate, but simply use the data given.
vi. .HVf. means the formation energy of vacancy and .HB. means binding energy (4 configurations will be used for .HB., so 4 words or 1 float are expected after .HB.). 
vii.    Current code will not likely to work if these keywords are spelt incorrectly.

Below are two examples of the $freq part in the input file:
Ex1:
$freq
type 5

v1 vac1 vac10-vac1 
v2  2
v3 vac3 vac4-vac3 
v4  5 
v0 vac0 vac00-vac0 

E1  vac1 vac10-vac1

E2 vac2   vac20-vac2 
E3   0.5
E4 vac4 vac4-vac3 
E0 vac0 vac00-vac0
HVf  0.5
#TTM can also have HVf  perfect   vac   ???
HB  perfect  sub  vac-sub  vac
$end


Ex2:
$freq

type hcp

HVf 0.44
HB -0.1
Ea 0.5
Eb 0.5 
Ec 0.5 
EX 0.5 
Eap 0.5
Ebp 0.5
Ecp 0.5
EXp 0.5
va 5 
vb 5 
vc 5
vX 3
vap 5
vbp 5
vcp 3
vXp 4

 
$end

**************************
Defect finder
**************************

The defect finder takes a POSCAR file and finds vacancies and interstitials.
The defect finder currently exists in a separate repository.
You may test it online at materialshub.org > Resources > Tools > Defect Finder


