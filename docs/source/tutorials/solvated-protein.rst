.. _tutorialSolvatedProtein:

Classical MD of a solvated protein
==================================

This tutorial takes the mini-protein Trp-cage (PDB `1L2Y
<https://www.rcsb.org/structure/1L2Y>`_) from the PDB file to a running PQ
simulation. It uses the AMBER ff14SB force field and TIP3P water:

1. OpenMM adds hydrogens, water and a chloride ion and minimises the system.
2. ParmEd writes an AMBER prmtop and inpcrd.
3. ``parm2pq`` turns them into PQ setup files.
4. PQ runs NVT and then NPT molecular dynamics.
5. PQAnalysis converts the energy file and computes a radial distribution
   function.

All files are in ``examples/parm2pq`` of the repository and are shown below.

.. admonition:: Before you start
   :class: warning

   Read the limits of ``parm2pq`` in :ref:`forceFieldImport` first. In short:
   PQ has no PME, so this tutorial uses a damped shifted-force (Wolf) Coulomb
   cutoff; CHARMM protein force fields (CMAP) are not supported; hydrogen
   mass repartitioning is lost.

What you need
-------------

- Python 3.12 or newer.
- A PQ build: follow the `PQ installation guide
  <https://molarverse.github.io/PQ/>`_. The commands below call it as ``PQ``.
- About one hour of one CPU core for the MD runs (PQ runs MM on one core).

Create an environment with PQAnalysis, ParmEd and OpenMM:

.. code-block:: console

   $ python -m venv pq-tutorial
   $ source pq-tutorial/bin/activate
   $ pip install "PQAnalysis[protein]"

``PQAnalysis[ff]`` installs ParmEd, which ``parm2pq`` needs. ``parm2pq`` is
new: if ``parm2pq --help`` fails, install the development version with
``pip install "PQAnalysis[ff] @ git+https://github.com/MolarVerse/PQAnalysis.git@dev"``.

Get the files
-------------

Copy the example folder and download the structure:

.. code-block:: console

   $ git clone --depth 1 https://github.com/MolarVerse/PQAnalysis.git
   $ cp -r PQAnalysis/examples/parm2pq trpcage
   $ cd trpcage
   $ curl -O https://files.rcsb.org/download/1L2Y.pdb

The folder holds the build script, three PQ and PQAnalysis input files and the
tleap files for the alternative route at the end.

Step 1: build the solvated system
---------------------------------

.. code-block:: console

   $ python build_system.py 1L2Y.pdb

The script does five things:

- It takes the first of the 38 NMR models and adds hydrogens for pH 7.
- It puts the protein in a cubic TIP3P box with 1 nm of water on each side
  and adds one Cl\ :sup:`-` to cancel the +1 charge of Trp-cage.
- It builds the system **without constraints and with flexible water**, so
  every bond keeps its force constant. ``parm2pq`` needs them; it adds the
  constraints for PQ itself.
- It minimises the energy with OpenMM, so PQ starts from a relaxed structure
  (see `Minimisation`_ for why not in PQ).
- It saves ``system.prmtop`` and ``system.inpcrd`` with ParmEd; the inpcrd
  holds the box.

``build_system.py``:

.. literalinclude:: ../../../examples/parm2pq/build_system.py
   :language: python

The script prints the returned preparation summary, including the output
paths, atom and water counts, box lengths in Angstrom, and minimized
OpenMM energy in kcal/mol. It will not replace an existing AMBER pair.

The minimised energy changes a little from run to run, because the
minimiser runs on several threads. It takes about 15 seconds.

Step 2: convert with parm2pq
----------------------------

.. code-block:: console

   $ parm2pq system.prmtop -c system.inpcrd -d md --shake-h

- ``-c`` gives the coordinates and the box.
- ``-d md`` writes the PQ files into the new folder ``md``.
- ``--shake-h`` turns every bond to hydrogen into a SHAKE constraint, so
  PQ can use a 2 fs timestep. The water O-H bonds are fixed, its H-O-H
  angle stays flexible.

``md`` now holds:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - File
     - Size
     - Content
   * - ``moldescriptor.dat``
     - 10 kB
     - three moltypes: the protein (304 atoms), water and chloride
   * - ``topology.top``
     - 107 kB
     - bonds, angles, dihedrals, impropers and 3164 SHAKE constraints
   * - ``parameter.prm``
     - 6 kB
     - force constants, 1-4 factors and Lennard-Jones pairs
   * - ``intra_nonbonded.dat``
     - 174 kB
     - excluded and scaled 1-4 pairs inside the protein and the water
   * - ``start.rst``
     - 467 kB
     - box, positions and moltypes of all 4826 atoms
   * - ``pq_input.in``
     - 0.4 kB
     - the input keys for these files, a starting point only

The system has a box, so ``pq_input.in`` contains ``long_range = wolf``
with ``wolf_param = 0.2``, as do the input files below; see
`Long-range Coulomb`_.

Step 3: the PQ input files
--------------------------

The run has two parts. A 10 ps NVT run heats the system and lets the
solvent relax at constant volume. A 20 ps NPT run then lets the box adjust
to 1 bar. Copy the inputs into ``md``:

.. code-block:: console

   $ cp nvt.in npt.in rdf.in md/
   $ cd md

``nvt.in``:

.. literalinclude:: ../../../examples/parm2pq/nvt.in
   :language: text

``npt.in``:

.. literalinclude:: ../../../examples/parm2pq/npt.in
   :language: text

What the keys do:

- ``jobtype = mm-md`` and ``force-field = on`` run classical MD with the
  bonded terms from ``topology.top`` and ``parameter.prm``.
- ``file_prefix`` names all output files (``nvt.en``, ``nvt.xyz``,
  ``nvt.rst`` and so on).
- ``rcoulomb = 12`` is the cutoff in Angstrom for Coulomb and
  Lennard-Jones. Each box length must be more than twice this value.
- ``cell-list = on`` finds the pairs within the cutoff with a cell list
  instead of checking all pairs.
- ``shake = on`` applies the constraints that ``--shake-h`` wrote.
- ``timestep = 2.0`` is in fs; it is safe because the fastest motions
  (bonds to hydrogen) are constrained.
- ``init_velocities = true`` draws start velocities at ``temp``; the
  start file from ``parm2pq`` has none. Without this key the run starts at
  0 K.
- ``thermostat = velocity_rescaling`` is the stochastic Bussi thermostat;
  it gives the canonical ensemble.
- ``manostat = stochastic_rescaling`` scales the box towards ``pressure``
  (bar) with the compressibility of water (the default).
  ``virial = molecular`` (also the default) computes the pressure from the
  molecular centres of mass.
- ``npt.in`` starts from ``nvt.rst``, the restart file the NVT run writes,
  so it keeps positions, velocities and the box.

Long-range Coulomb
~~~~~~~~~~~~~~~~~~

PQ has no PME. With ``long_range = none`` (the PQ default) it cuts the
Coulomb interaction at ``rcoulomb`` and only shifts it to zero there. In a box of
water with ions, the interactions beyond 12 Angstrom are then dropped
without damping, which distorts the structure of the solvent and the ions.
``long_range = wolf`` damps the interaction with an error function before
the cutoff (damped shifted force). With ``wolf_param = 0.2`` per Angstrom
and a 12 Angstrom cutoff, energies and forces are close to an Ewald sum
(Fennell and Gezelter, J. Chem. Phys. 124, 234104, 2006).
``long_range = reaction-field`` with ``rf_epsilon = 78.5`` is a cheaper
alternative for water.

Minimisation
~~~~~~~~~~~~

The system is minimised in OpenMM before ``parm2pq``, not with PQ's
``mm-opt``. ``mm-opt`` ignores SHAKE, so after ``--shake-h`` the hydrogen
atoms have no bond and fly off within a few steps. Without ``--shake-h`` it
is stable, but it needs ``learning-rate-strategy = constant`` (the default
strategy needs ``learning-rate-decay``, which has no default) and it ends
with ``Error: Optimizer did not converge`` once ``nstep`` is reached. Always
minimise before ``parm2pq``, with OpenMM or ``sander``.

Step 4: run PQ
--------------

.. code-block:: console

   $ PQ nvt.in
   $ PQ npt.in

Each run writes its files with the prefix from ``file_prefix``: energies
(``.en``, column names in ``.info``), trajectory (``.xyz``), velocities
(``.vel``), box (``.box``), the restart file (``.rst``) and the log. A PQ
build with Python support may print ``Could not find platform dependent
libraries <exec_prefix>`` at the start; it does not affect MM runs.

Check that the runs ended cleanly: ``nvt.log`` and ``npt.log`` end with
``PQ ended normally``.

What to expect, from a test run on one core of an Intel Core i7-11700
(Arch Linux, btrfs; the NPT numbers are from its first 5 ps):

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * -
     - NVT, 10 ps
     - NPT, first 5 ps
   * - time per step
     - 0.29 s (24 min in total)
     - 0.19 s
   * - temperature, second half
     - 299.5 +/- 2.7 K
     - 301.0 +/- 2.0 K
   * - pressure, second half
     - -1110 +/- 230 bar
     - -150 +/- 220 bar
   * - box length
     - 37.50 Angstrom
     - 37.50 to 36.55 Angstrom
   * - density of the whole system
     - 0.92 g cm\ :sup:`-3`
     - 0.93 to 1.00 g cm\ :sup:`-3`

The water box from OpenMM starts below the density of liquid water, so the
NVT pressure is strongly negative. The NPT run shrinks the box and the
pressure moves towards 1 bar. This is why the NPT step is needed before any
production run. The total momentum stays below 10\ :sup:`-12` amu Angstrom
fs\ :sup:`-1`.

Step 5: analyse
---------------

Convert the NPT energy file into an XVG table; ``--timestep 2`` turns the
step column into ps:

.. code-block:: console

   $ pqanalysis convert npt.en -o npt.xvg --timestep 2

``npt.xvg`` holds every column of the energy file (temperature, pressure,
the energy terms, momentum and time per step) and opens in xmgrace. By
default it plots four energies; choose the columns with ``--y``:

.. code-block:: console

   $ pqanalysis convert npt.en -o temperature.xvg --timestep 2 --y TEMPERATURE --y PRESSURE

The box lengths of every output frame are in ``npt.box``.

Compute the oxygen-oxygen RDF of the water:

.. literalinclude:: ../../../examples/parm2pq/rdf.in
   :language: text

.. code-block:: console

   $ pqanalysis rdf rdf.in

``304..3..4822`` selects every third atom from 304 to 4822, the water
oxygens: the protein has atoms 0 to 303, then come the 1507 waters (O, H,
H) and the chloride ion. Indices start at 0; see :doc:`../data/selections`.
If ``build_system.py`` printed another number of waters :math:`N`, the last
oxygen is :math:`301 + 3N`.

The result is in ``rdf_OO.dat`` (columns described in
:ref:`analysis-output-rdf`). In the test run the first peak is at
2.78 Angstrom with :math:`g(r) \approx 2.9` and the first minimum is near
3.5 Angstrom, as expected for TIP3P water.

Alternative: build the system with tleap
----------------------------------------

AmberTools gives the same kind of prmtop and inpcrd. It needs a PDB with
AMBER residue names and the first model only:

.. code-block:: console

   $ pdb4amber -i 1L2Y.pdb -o 1L2Y_amber.pdb --model 1 --nohyd
   $ tleap -f tleap.in
   $ sander -O -i min.in -p system.prmtop -c system.inpcrd -r system_min.rst7
   $ parm2pq system.prmtop -c system_min.rst7 -d md --shake-h

``tleap.in`` loads ff14SB and TIP3P, solvates with 10 Angstrom of water in a
near-cubic box and adds the chloride ion:

.. literalinclude:: ../../../examples/parm2pq/tleap.in
   :language: text

``min.in`` minimises the system with ``sander``:

.. literalinclude:: ../../../examples/parm2pq/min.in
   :language: text

From here, continue with step 3, with one change. tleap's TIP3P water has
an H-H bond, so it is rigid: ``parm2pq`` writes ``mshake.dat`` with one
reference water and keeps the waters rigid with M-SHAKE, while the protein
hydrogens stay in SHAKE. Replace ``shake = on`` in ``nvt.in`` and ``npt.in``
with the two lines from ``pq_input.in``:

.. code-block:: text

   shake = mshake;
   mshake_file = mshake.dat;

The tleap box is larger (about 49 Angstrom, 2769 waters, 8612
atoms), so each step takes longer. tleap puts the chloride ion before the
water, so the water oxygens in ``rdf.in`` are ``305..3..8609``.

Next steps
----------

- Make the NPT run longer for production (change ``nstep``) and analyse
  the protein with :doc:`../analyses/msd` or more RDFs.
- For other systems, check the supported terms and the troubleshooting list
  in :ref:`forceFieldImport`.
