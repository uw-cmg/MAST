*****************************
The Scaling section
*****************************

The ``$scaling`` section supports finite-size scaling through different supercell sizes. 

================================
Finite-size scaling
================================
Finite size scaling is supported with a special "scaling" subsection.

Defect positions will be automatically scaled.

    * For example, ``0.25 0.0 0.0`` in the original supercell would become ``0.125 0.0 0.0`` in a 2x1x1 cell. 

Special notes:

*  :doc:`3_1_2_ingredients` should include an "inducescaling" ingredient with a ``mast_run_method`` of ``run_scale``

*  :doc:`3_1_3_recipe` should include ``inducescaling_<S>`` and ``defect_<S>`` ingredients.

    *  The "<S>" tags will correspond to the scaling sizes and labels.

For each scaling size, create a subsection with ``begin <labelname>``.

Within the subsection, include keywords:

    * **mast_size** with a scaling matrix of integers ``[M, N, P]`` or ``[M1 M2 M3, N1 N2 N3, P1 P2 P3]``
    
    * **mast_kpoints** with a Kpoint mesh in the form ``QxRxS``, followed by a 

    * Kpoint mesh type, M for Monkhorst-Pack and G for Gamma-point centered

    * (Optional) Kpoint mesh shift, in floats, e.g. ``0.1 0.2 0.3`` 
    
    * **WILL OTHER KEYWORDS NOW WORK AS WELL, FOR OVERRIDING?**

Example::
 
    $scaling
    begin 2xhigh
    mast_size [2, 2, 2]
    mast_kpoints 4x4x4 M
    end
    begin 4xhigh
    mast_size [4, 4, 4]
    mast_kpoints 2x2x2 M
    end
    $end

**NEED TO VERIFY MADELUNG POTENTIAL CHANGES**

In order to figure out which scaling sizes to use for finite-size scaling, MAST includes a Madelung potential utility.

This utility generates a distribution of cell sizes for best scaling, according to the method in::

    Hine, N. D. M., Frensch, K., Foulkes, W. M. C. & Finnis, M. W. Supercell size scaling of density functional theory formation energies of charged defects. Physical Review B 79, 13, doi:10.1103/PhysRevB.79.024112 (2009).

Run this utility as follows in order to generate a cut-and-paste for the scaling section. ::

    mast_finite_size_scaling_sizes perfDir defDir minDefDist maxNumAtoms numStructAsked

* **perfDir**: perfect primordial (small) cell directory, which should already have run and include VASP CONTCAR, OSZICAR, etc. files.

* **defDir**: defected primordial cell directory, which should already have run and include VASP CONTCAR, OSZICAR, etc. files.

* **minDefDist** (default 3): minimum defect-defect distance between periodic images, in Angstroms.

* **maxNumAtoms** (default 600): maximum number of atoms for supercell size evaluations

* **numStructAsked** (default 5): number of structures to return in the distribution 

* Note that you will have to manually adjust the kpoint mesh in your cut-and-paste.


.. raw:: html

    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-54660326-1', 'auto');
      ga('send', 'pageview');

    </script>

