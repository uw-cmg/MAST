Input files for GA+STEM simulation 
-------------------------------------

Au_u3.eam: EAM potential file

PSF.txt: pre-calculated PSF, 976 x 976 pixels  

STEM_ref: file including atomic coordinates of a Au309 Ino-Decahedron structure

STEM_run.py: set up parameters for simualtion

submit.sh: submit script


Output files for GA+STEM simulation
-----------------------------------

Output-rank0.txt: main output file, include all information through evolution 

Output.txt: from 2 to 5 cols: generation #; fitness function; energy per atom in eV; alpha*chi2

Output-rank0/Summary-Output.txt: summary of fitness evolution and computing cost 


Visulization
------------

image.py: readin a strucutre (xyz format) and precalculated PSF, simulate STEM image using convolution method.  
