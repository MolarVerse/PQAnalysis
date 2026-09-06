Example data
============

The repository ships a checkout-only tutorial fixture under
``examples/water``. It is an isolated water molecule in a 10 Å cubic cell
(25 frames) plus the bundled H₂O Hessian. Use it to copy-paste the Python and
CLI recipes. Do not treat the numerical results as bulk-liquid RDF, diffusion
or infrared data.

Get the files
-------------

Clone the repository (the fixture is not installed by pip):

.. code-block:: console

   $ git clone https://github.com/MolarVerse/PQAnalysis.git
   $ cd PQAnalysis/examples/water

After this pull request lands, the same folder is on the default branch.

What is in the folder
---------------------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - File
     - Role
   * - ``trajectory.xyz``
     - 25 PQ XYZ frames of one water molecule (O, H, H)
   * - ``trajectory.vel``
     - Matching velocity frames for VACF and momentum
   * - ``restart.rst``
     - Restart of the first frame (topology / residue ids)
   * - ``moldescriptor.dat``
     - H₂O residue template with partial charges
   * - ``structure.rst`` / ``hessian.dat``
     - Isolated-molecule Hessian fixture for vibrations
   * - ``rdf.in``, ``msd.in``, ``vacf.in``, ``vibrations.in``
     - Input files sized to this 25-frame trajectory

``msd.in`` and ``vacf.in`` use ``window = 8``. Production trajectories need
much larger windows; those keys are documented on the method pages.

Run the first RDF
-----------------

.. code-block:: console

   $ pqanalysis rdf rdf.in

The table begins (values will match bit-for-bit on this fixture):

.. code-block:: text

   # PQAnalysis: Radial distribution function
   # FIELDS r_i g_r_i N_r_i g_r_i_dV_i H_i_minus_E_i
   0.25 0.0 ...
   0.75 261.9 ...

The large :math:`g(r)` near 0.75 Å is the intramolecular O–H peak of a single
molecule in a large box, not a liquid first shell. See
:ref:`analysis-output-rdf` for column definitions.

The other commands in this folder are:

.. code-block:: console

   $ pqanalysis msd msd.in
   $ pqanalysis vacf vacf.in
   $ pqanalysis vibrations vibrations.in
   $ pqanalysis check_momentum trajectory.vel \
       --selection all --output momentum.dat

Python recipes that load these files are on :doc:`python-api`. A continuous
load → RDF → MSD → figure script is on :doc:`workflow`.
