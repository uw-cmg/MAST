###################################
The Defects section
###################################
The ``$defects`` section specifies defects by:

*  defect type:
    *  vacancy
    *  interstitial
    *  substitution or antisite

*  defect coordinates

*  defect element symbol
    *  Note that if an ``elementmap`` subsection is given in :doc:`3_1_1_structure`, then the mapped designations ``X1``, ``X2``, and so on can be given instead of an element symbol.

**ATTENTION:** 

*  Elements in the initial structure, given in :doc:`3_1_1_structure`, will appear in order as entered, by posfile keyword or through the coordinates and/or elementmap subsections.

*  However, once a defect is formed, structures are RESORTED by element ELECTRONEGATIVITY.  Therefore, if you are using substitutional defects or non-self-interstitials, you may find that later element-specific keywords (mast_setmagmom, LDAUU, LDAUJ) may be OUT OF ORDER FOR YOUR DEFECTED STRUCTURE.

*  Please check your files carefully! You may want a separate input file for each chemical system (possibly created through looping (see :doc:`3_0_inputfile`) in order to synchronize the elements and element-specific keywords.

The ``coord_type`` keyword specifies fractional or cartesian coordinates for the defects.

The ``threshold`` keyword specifies the absolute threshold for finding the defect coordinate, since relaxation of the perfect structure may result in changed coordinates.

Example ``$defects`` section::

    $defects

    coord_type fractional
    threshold 1e-4

    vacancy 0 0 0 Mg
    vacancy 0.5 0.5 0.5 Mg
    interstitial 0.25 0.25 0 Mg
    interstitial 0.25 0.75 0 Mg
    
    $end

The above section specifies 4 point defects (2 vacancies and 2 interstitials) to be applied separately and independently to the structure. When combined with the correct recipe in :doc:`3_1_3_recipe`, four separate ingredients, each containing one of the defects above, will be created.

Multiple point defects can be also grouped together as a combined defect within a ``<defect label>`` subsection such as::

    $defects
    
    coord_type fractional
    threshold 1e-4
    
    begin doublevac
    vacancy 0.0 0.0 0.0 Mg
    vacancy 0.5 0.5 0.5 Mg
    end
    
    interstitial 0.25 0.25 0 Mg
    interstitial 0.25 0.75 0 Mg
    
    $end

In this case, there will be three separate defect ingredients: one ingredient with two vacancies together (where the defect group is labeled ``doublevac``), one interstitial, and another interstitial.


=====================
Charges for defects
=====================
Charges can be specified as ``charge=0,10``, where a comma denotes the lower and upper ranges for the charges.

Let's say we want a Mg vacancy with charges from 0 to 3 (0, 1, 2, and 3)::

    vacancy 0 0 0 Mg charge=0,3

Let's say we want a dual Mg vacancy with a charge from 0 to 3 and labeled as Vac@Mg-Vac@Mg::

    begin Vac@Mg-Vac@Mg
    vacancy 0.0 0.0 0.0 Mg
    vacancy 0.5 0.5 0.5 Mg
    charge=0,3
    end

For a single defect, charges and labels can be given at the same time:

Let's say we have a Mg vacancy with charges between 0 and 3, and we wish to label it as Vac@Mg::

    vacancy 0.0 0.0 0.0 Mg charge=0,3 label=Vac@Mg

The charge and label keywords are interchangeable, i.e. we could also have typed::

    vacancy 0 0 0 Mg label=Vac@Mg charge=0,3

If you use charges in the defects section like this, then you must use a tagged ``defect_<N>_<Q>`` type recipe in :doc:`3_1_3_recipe`.

=====================
Phonons for defects
=====================

Phonon calculations are described by a *phonon center site* coordinate and a *phonon center radius* in Angstroms. Atoms within the sphere specified by these two values will be included in phonon calculations.

For VASP, this inclusion takes the form of selective dynamics T T T for the atoms within the sphere, and F F F otherwise, in a phonon calculation (IBRION = 5, 6, 7, 8)

If the phonon center radius is 0, only the atom found at the phonon center site point will be considered.

To use phonons in the defects section, use the subsection keyword ``phonon`` followed by:

* A label for the phonon

* The fractional coordinates for the phonon center site

*  A float value for the phonon center radius

*  An optional float value for the tolerance-matching threshold for matching the phonon center site (if this last value is not specified, 0.1 is used). 

Multiple separate phonon calculations may be obtained for each defect, for example::

    begin int1
    interstitial 0.25 0.25 0.25 X2
    phonon host3 0.3 0.3 0.4 2.5 0.01
    phonon solute 0.1 0.1 0.2 0.5
    end

In the example above, *host3* is the label for the phonon calculation where (0.3, 0.3, 0.4) is the coordinate for the phonon center site, and 2.5 Angstroms is the radius for the sphere inside which to consider atoms for the phonon calculation. Points within 0.01 of fractional coordinates will be considered for matching the phonon center site. 

In the example above, *solute* is the label for the phonon calculation bounded within a 0.5 Angstrom radius centered at (0.1, 0.1, 0.2) in fractional coordinates. As no threshold value was given, points within 0.1 (default) of fractional coordinates will be considered for matching the phonon center site.


The recipe template file for phonons may include either the explicit phonon labels and other labels, or <S>, <N>, <Q>, <P>. See :doc:`3_1_3_recipe`.

Because phonons are cycled with the defects, a new parent loop must be provided for the phonons, for example::

    {begin}
    defect_<N>_<Q>_stat (static)
        phonon_<N>_<Q>_<P> (phonon)
            phonon_<N>_<Q>_<P>_parse (phononparse)
    {end}

