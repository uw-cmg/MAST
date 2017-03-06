#########################################
Particle Trajectory Diffusion Analysis
#########################################
Author: Leland Barnard
Acknowledgments to: Amy Bengtson, Saumitra Saha

Computes mean squared displacements and tracer diffusion coefficients from particle position data as a function of time.

`Version 1.1 - originally published on 28 Mar 2014 <https://nanohub.org/resources/22854>`_

* You must download and unzip the MAST tar.gz file from https://github.com/uw-cmg/MAST/releases in order to access the source code, which is in MAST/standalone/diffanalyzer, and will likely need to be recompiled for your machine. 

This tool takes as input particle position data from methods such as molecular dynamics or kinetic Monte Carlo and computes the mean squared displacement for all particles as a function of time. 
For a system with multiple types of particles, the mean squared displacement is computed for each particle type. 
The tracer diffusion coefficient is then calculated from the slope of the mean squared displacement vs time curve.

The tool is based on *The Working Man's Guide to Obtaining Self Diffusion Coefficients from Molecular Dynamics Simulations* by Professor David Keffer from UT Knoxville.

====================
Input file format:
====================

This tool reads in atomic position data in the VASP XDATCAR format. This file format begins with the following set of lines::

    Name
    C1 C2 C3 ...
    N1 N2 N3 ...
    Direct

* The first line is a name or description of the file. It is not read by the tool. 
* The second line are the names of the components in the system. These will be element names in the case of an atomic simulation. 
* The third line are the number of particles of each component in the system. 
* The final line is a VASP generated line that specifies direct atomic coordinates. 
* Following these 4 lines, the file must have 1 blank line, and then the particle position data begins on line 6. Particle positions must be in fractional or direct coordinates, and a single line must separate the blocks of particle positions at each time step throughout the file.

===================================================
Calculation of error on the diffusion coefficient:
===================================================

The error bars on the mean squared displacements represent a single standard deviation in the measurements of the squared displacements over all time origins.

The error in the diffusion coefficient represents the standard error in the slope of the weighted least squares fit to the mean squared displacement, using the variance in the squared displacements as the error weight.

============
References
============
"The Working Man's Guide to Obtaining Self Diffusion Coefficients from Molecular Dynamics Simulations" by Professor David Keffer from UT Knoxville, which may be found here: `<http://www.cs.unc.edu/Research/nbody/pubs/external/Keffer/selfD.pdf>`_

===============
Cite this work
===============
Researchers should cite this work as follows::

    Leland Barnard, Amy Bengtson, Saumitra Saha, Weixi Ma, Navjeet Sandhu, Henry Wu, and Dane Morgan, "Particle Trajectory Diffusion Analysis Tool," https://nanohub.org/resources/diffanalyzer (2015).
