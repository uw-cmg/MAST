================================================
Effective Grain Boundary Diffusivity Calculator
================================================
Author: Jie Deng

Calculates the effective diffusivity in a grain boundary network with two types of randomly distributed grain boundaries.

`Version 1.3 - originally published on 21 Feb 2014 <https://nanohub.org/resources/22858>`_

* You must download and unzip the MAST tar.gz file from https://github.com/uw-cmg/MAST/releases in order to access the source code, which is in MAST/standalone/gbdiff, and will likely need to be recompiled for your machine. 

This tool calculates the effective diffusivity in a grain boundary network represented by a three-dimensional Voronoi diagram. 
Two types of grain boundaries with different diffusivities are randomly distributed in the domain. 
The effective diffusivity is calculated using the mean squared displacement method, where periodic boundary conditions are applied in all directions. 
Users are free to choose the fraction of each grain boundary type as well as the activation energy and pre-factor for each grain boundary diffusivity.

=================
Cite this work
=================
Researchers should cite this work as follows::

    Jie Deng, Navjeet Sandhu, Henry Wu, Dane Morgan, and Izabela Szlufarska, "Effective Grain Boundary Diffusivity Calculator," https://nanohub.org/resources/gbdiffusion (2015).

.. raw:: html

    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-54660326-1', 'auto');
      ga('send', 'pageview');

    </script>

