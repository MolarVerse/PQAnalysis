.. _references:

References
==========

The entries below are the primary sources for the estimators PQAnalysis
implements. Each analysis page cites the entries that define the quantity it
reports, so a documented estimator can be checked against its original
definition rather than against this implementation alone.

Molecular simulation and liquid-state theory
--------------------------------------------

.. [Allen2017] Allen, M. P.; Tildesley, D. J. *Computer Simulation of
   Liquids*, 2nd ed.; Oxford University Press: Oxford, 2017.
   ISBN 978-0-19-880319-5.
   `doi:10.1093/oso/9780198803195.001.0001
   <https://doi.org/10.1093/oso/9780198803195.001.0001>`__

.. [Frenkel2002] Frenkel, D.; Smit, B. *Understanding Molecular Simulation:
   From Algorithms to Applications*, 2nd ed.; Academic Press: San Diego,
   2002. ISBN 978-0-12-267351-1.
   `doi:10.1016/B978-0-12-267351-1.X5000-7
   <https://doi.org/10.1016/B978-0-12-267351-1.X5000-7>`__

.. [Hansen2013] Hansen, J.-P.; McDonald, I. R. *Theory of Simple Liquids:
   with Applications to Soft Matter*, 4th ed.; Academic Press: Oxford, 2013.
   ISBN 978-0-12-387032-2.
   `doi:10.1016/C2010-0-66723-X <https://doi.org/10.1016/C2010-0-66723-X>`__

Transport coefficients and time correlation functions
------------------------------------------------------

.. [Einstein1905] Einstein, A. Über die von der molekularkinetischen Theorie
   der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten suspendierten
   Teilchen. *Annalen der Physik* **1905**, *322* (8), 549-560.
   `doi:10.1002/andp.19053220806
   <https://doi.org/10.1002/andp.19053220806>`__

.. [Green1954] Green, M. S. Markoff Random Processes and the Statistical
   Mechanics of Time-Dependent Phenomena. II. Irreversible Processes in
   Fluids. *The Journal of Chemical Physics* **1954**, *22* (3), 398-413.
   `doi:10.1063/1.1740082 <https://doi.org/10.1063/1.1740082>`__

.. [Kubo1957] Kubo, R. Statistical-Mechanical Theory of Irreversible
   Processes. I. General Theory and Simple Applications to Magnetic and
   Conduction Problems. *Journal of the Physical Society of Japan* **1957**,
   *12* (6), 570-586.
   `doi:10.1143/JPSJ.12.570 <https://doi.org/10.1143/JPSJ.12.570>`__

.. [Rahman1964] Rahman, A. Correlations in the Motion of Atoms in Liquid
   Argon. *Physical Review* **1964**, *136* (2A), A405-A411.
   `doi:10.1103/PhysRev.136.A405
   <https://doi.org/10.1103/PhysRev.136.A405>`__

Spectra from time correlation functions
----------------------------------------

.. [Wiener1930] Wiener, N. Generalized harmonic analysis. *Acta Mathematica*
   **1930**, *55*, 117-258.
   `doi:10.1007/BF02546511 <https://doi.org/10.1007/BF02546511>`__

.. [Khintchine1934] Khintchine, A. Korrelationstheorie der stationären
   stochastischen Prozesse. *Mathematische Annalen* **1934**, *109* (1),
   604-615.
   `doi:10.1007/BF01449156 <https://doi.org/10.1007/BF01449156>`__

.. [Dickey1969] Dickey, J. M.; Paskin, A. Computer Simulation of the Lattice
   Dynamics of Solids. *Physical Review* **1969**, *188* (3), 1407-1418.
   `doi:10.1103/PhysRev.188.1407
   <https://doi.org/10.1103/PhysRev.188.1407>`__

.. [Thomas2013] Thomas, M.; Brehm, M.; Fligg, R.; Vöhringer, P.; Kirchner, B.
   Computing vibrational spectra from ab initio molecular dynamics.
   *Physical Chemistry Chemical Physics* **2013**, *15* (18), 6608-6622.
   `doi:10.1039/C3CP44302G <https://doi.org/10.1039/C3CP44302G>`__

.. [Harris1978] Harris, F. J. On the use of windows for harmonic analysis
   with the discrete Fourier transform. *Proceedings of the IEEE* **1978**,
   *66* (1), 51-83.
   `doi:10.1109/PROC.1978.10837
   <https://doi.org/10.1109/PROC.1978.10837>`__

Normal modes and infrared intensities
--------------------------------------

.. [Wilson1955] Wilson, E. B., Jr.; Decius, J. C.; Cross, P. C. *Molecular
   Vibrations: The Theory of Infrared and Raman Vibrational Spectra*;
   McGraw-Hill: New York, 1955. Reprinted by Dover: New York, 1980.
   ISBN 978-0-486-63941-3.

.. [Person1974] Person, W. B.; Newton, J. H. Dipole moment derivatives and
   infrared intensities. I. Polar tensors. *The Journal of Chemical Physics*
   **1974**, *61* (3), 1040-1049.
   `doi:10.1063/1.1681972 <https://doi.org/10.1063/1.1681972>`__

.. _software-provenance:

Software provenance
-------------------

Several PQAnalysis kernels reproduce the arithmetic of older programs bit for
bit so that historical results remain reproducible. That provenance is
recorded here as software attribution, not as literature.

.. [thhTools] The ``thh_tools`` collection of legacy analysis programs: the
   ``RDF`` C code, the ``Diffcalc`` mean-square-displacement code, the
   ``FreqCalc`` and ``Fluxfreqcalc`` correlation codes, the ``ftvac``
   (``ft.f``) Fourier transformation code and ``thh_momentum/equipartition.jl``.
   The collection is not publicly distributed and has no release version or
   DOI. PQAnalysis reproduces its operation order only on the compatibility
   paths described on the individual analysis pages, and pins that behavior
   with the parity fixtures described in :doc:`developerGuide/validation`.

.. _citing-pqanalysis:

Citing PQAnalysis
-----------------

PQAnalysis has no journal article, no ``CITATION.cff`` file and no minted DOI.
Do not cite a DOI for it, because none exists.

Cite the software by repository and by the exact release you ran:

.. code-block:: text

   PQAnalysis (version <release>), MolarVerse.
   https://github.com/MolarVerse/PQAnalysis

Report the version string from ``pqanalysis --version`` or from
``PQAnalysis.__version__``. Reported analyses should also name the
simulation engine that produced the input data, normally
`PQ <https://github.com/MolarVerse/PQ>`__, and the method-specific settings
that change the reported numbers, such as bin width, correlation window,
apodization function or Hessian unit.

If a DOI is minted for a future release, cite that DOI instead of this
section.
