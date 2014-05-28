.. _8_0_1_diffusionconnectivity:

***************************
Diffusion Connectivity
***************************
A new section is introduced in the input file::

    $site
    int1
    0.5 0.5 0.5
    0.5 0 0
    0 0.5 0
    0 0 0.5
    int2
    0.25 0.25 0.25
    0.75 0.25 0.25
    0.25 0.75 0.25
    0.25 0.25 0.75
    0.75 0.75 0.75
    0.25 0.75 0.75
    0.75 0.25 0.75
    0.75 0.75 0.25 
    $end

In this example, there are two types of local minimum (interstitial site) int1 and int2. The geometrically equivalent site coordinates are listed for each type.

The create_paths.py code first parses the perfect lattice with the defect sites, then finds the nth neighbor of possible pairs of both same and different site types, and detect if the path candidate crosses over the host lattice site or another defect site and delete it. 
If the lattice user provides is too small and not all the neighbors (up to nth) are found, the code will double, triple, etc. the size until all required neighbor pairs are found.

The NEBcheck.py code will generate the MAST-input style defect structure for the possible pairs found in the create_paths.py and write a new input file that calls MAST to generate NEB and phonon folders to check if these paths are appropriate. Currently the code ends at calling MAST and does not yet manage to handle the NEB and phonon results analysis.

The usage is: python NEBcheck.py -i <input> -n <up-to-nth-neighbor> 
