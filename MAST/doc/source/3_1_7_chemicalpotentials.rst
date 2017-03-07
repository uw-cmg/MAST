###################################
The Chemical Potentials section
###################################

The ``$chemical_potentials`` section lists chemical potentials, used for defect formation energy calculations using the defect formation energy tool.

Currently, chemical potentials must be set ahead of time. Each chemical potential subsection may be labeled descriptively. ::

    $chemical_potentials
    
    begin Ga rich
    Ga -3.6080
    As -6.0383
    Bi -4.5650
    end
    
    begin As rich
    Ga -4.2543
    As -5.3920
    Bi -4.5650
    end
    
    $end


.. raw:: html

    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-54660326-1', 'auto');
      ga('send', 'pageview');

    </script>

