Analyses
========

PQAnalysis covers structural organization, translational dynamics,
time-correlation spectra, molecular normal modes and conservation diagnostics.
Choose the observable from the physical question and from the data recorded by
the simulation.

Choose by input data
--------------------

.. list-table:: Analysis inputs and primary observables
   :class: pq-record-table pq-observable-table
   :header-rows: 1
   :widths: 24 32 44

   * - Analysis
     - Required physical data
     - Primary observable
   * - RDF
     - Positions and periodic cell
     - :math:`g_{AB}(r)` and cumulative coordination
   * - MSD
     - Positions and periodic cell
     - :math:`\langle |\mathbf{r}(t)-\mathbf{r}(0)|^2\rangle`
   * - VACF
     - Velocities and frame time step
     - Normalized :math:`C_v(t)` and its spectrum
   * - Vibrations
     - Structure, masses and Cartesian Hessian
     - Normal-mode wavenumbers and force constants
   * - Momentum
     - Velocities and atomic masses
     - :math:`|\sum_i m_i\mathbf{v}_i|` per frame

The method pages define each estimator, its assumptions and its interpretation
limits. File columns and units are specified once in
:ref:`analysisOutputFiles`. Programmatic entry points are listed in the
:doc:`../reference/functions`.

.. toctree::
   :hidden:
   :maxdepth: 1

   rdf
   msd
   vacf
   vibrations
   momentum
