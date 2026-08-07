Analyses
========

PQAnalysis covers structural organization, translational dynamics,
time-correlation spectra, molecular normal modes and conservation diagnostics.
Choose the observable from the physical question and from the data recorded by
the simulation.

.. grid:: 1 2 3 3
   :gutter: 2

   .. grid-item-card:: Radial distribution
      :link: rdf
      :link-type: doc

      Pair structure, preferred separations and coordination numbers.

   .. grid-item-card:: Mean square displacement
      :link: msd
      :link-type: doc

      Translational motion and Einstein-relation diffusion estimates.

   .. grid-item-card:: VACF and spectra
      :link: vacf
      :link-type: doc

      Velocity or charge-flux correlation and frequency-domain spectra.

   .. grid-item-card:: Vibrational analysis
      :link: vibrations
      :link-type: doc

      Hessian normal modes, wavenumbers, force constants and IR intensities.

   .. grid-item-card:: Total momentum
      :link: momentum
      :link-type: doc

      Frame-resolved linear momentum and center-of-mass drift diagnostics.

   .. grid-item-card:: Output schemas
      :link: ../userGuide/analysisOutputFiles
      :link-type: doc

      Exact columns, units, normalizations and format conversion behavior.

Choose by input data
--------------------

.. list-table:: Analysis inputs and primary observables
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

.. toctree::
   :hidden:
   :maxdepth: 1

   rdf
   msd
   vacf
   vibrations
   momentum
