.. _forceFieldImport:

Force-Field Import (parm2pq)
============================

``parm2pq`` converts a parametrized system from AMBER, CHARMM or GROMACS into
the setup files of a PQ ``mm-md`` run. It reads the input with
`ParmEd <https://parmed.github.io/ParmEd/>`_, an optional dependency:

.. code-block:: console

   $ pip install PQAnalysis[ff]

For a full run from a PDB file to a PQ simulation, see the tutorial
:doc:`../tutorials/solvated-protein`.

.. admonition:: Before you start
   :class: warning

   - **Long-range Coulomb.** PQ has no PME or Ewald sum. ``parm2pq`` writes
     ``long_range = none`` into ``pq_input.in``: a plain shifted-force
     cutoff, the form used in the validation below. In a periodic box with
     water and ions it drops all charge-charge interactions beyond the
     cutoff without damping, which distorts the structure of ions and water. For periodic
     charged systems use ``long_range = wolf`` with ``wolf_param = 0.2``
     (damped shifted force, a good match for a 12 Angstrom cutoff), or
     ``long_range = reaction-field`` with ``rf_epsilon = 78.5`` for water.
   - **CHARMM proteins are not supported yet.** CHARMM protein force fields
     need CMAP corrections, which PQ cannot evaluate, so ``parm2pq`` stops
     with an error. CHARMM files without CMAP, Urey-Bradley terms and
     harmonic impropers convert.
   - **Hydrogen mass repartitioning is lost.** PQ takes every mass from the
     element, so a repartitioned topology runs with normal masses. Use
     ``--shake-h`` and a 2 fs timestep, not 4 fs.
   - **Only AMBER is validated against OpenMM.** AMBER energies and forces
     match OpenMM (see `Validation`_). GROMACS and CHARMM input is only
     tested by round trips through ParmEd files.

Usage
-----

.. code-block:: console

   $ parm2pq system.prmtop -c system.inpcrd -d pq_setup
   $ parm2pq topol.top -c conf.gro -d pq_setup --shake-h
   $ parm2pq system.psf -c system.crd -p top.rtf par.prm -d pq_setup

The output directory receives:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Content
   * - ``moldescriptor.dat``
     - one moltype per unique molecule: element, atom type index, charge
       and global van der Waals index of each atom
   * - ``topology.top``
     - bonds, angles, dihedrals, impropers and (``--shake-h``) SHAKE
       constraints by global atom index
   * - ``parameter.prm``
     - the 1-4 scaling factors, the bonded parameters and one
       ``NONCOULOMBICS LJ`` line per pair of van der Waals types
   * - ``intra_nonbonded.dat``
     - the non-bonded pairs inside each moltype; negative partners are
       scaled 1-4 pairs, excluded pairs are left out
   * - ``start.rst``
     - box, element, moltype, positions and, if the coordinate file has
       them, velocities
   * - ``pq_input.in``
     - the input keys that load these files (``force-field = on``,
       ``rcoulomb``, ``long_range = none`` and the file names); add the run
       keys such as ``nstep``, ``timestep`` and the thermostat
   * - ``mshake.dat``
     - only with ``--shake-h`` and rigid water: the reference geometry of
       each rigid water moltype for M-SHAKE

Existing files are only replaced with ``--overwrite``. A structure without a
box gets a cubic box of its extent plus two Coulomb cutoffs, so no atom sees a
periodic image (vacuum).

Mapping
-------

ParmEd stores every force field in AMBER units and conventions.

- **Molecules.** Every connected group of bonded atoms is one molecule;
  protein chains joined by a bond (for example a disulfide) are split by
  chain. Identical molecules share a moltype. PQ takes the mass of each atom
  from its element, so hydrogen mass repartitioning is not kept.
- **Bonds and angles.** PQ uses :math:`\frac{1}{2} k`, AMBER :math:`k`, so
  the force constant is doubled.
- **Proper and periodic improper torsions.** One parameter line per Fourier
  term. PQ evaluates :math:`k (1 + \cos(n \varphi + \delta))`, AMBER
  :math:`k (1 + \cos(n \varphi - \delta))` with the same sign of
  :math:`\varphi`, so the phase changes sign.
- **Lennard-Jones.** The pair parameters :math:`c_6 = -B` and
  :math:`c_{12} = A` come from the AMBER tables or from the combining rule
  and NBFIX entries, so pair-specific parameters are kept.
- **Exclusions and 1-4 pairs.** Atoms one or two bonds apart and the ends
  of proper torsions are excluded, the ends of torsions that keep their 1-4
  interaction (or the GROMACS ``[ pairs ]``) are scaled. Inside a molecule
  this goes into the intra non-bonded file. Between molecules, the first
  bond, angle or proper torsion that connects such a pair carries the linker
  mark ``*``, so PQ removes (or scales) the pair it computes as
  intermolecular.
- **Constraints.** ``--shake-h`` turns every bond to hydrogen into a
  SHAKE constraint. A water with an H-H bond (as in tleap's TIP3P) is
  rigid: its two O-H bonds, the H-H bond and the angle are replaced by one
  M-SHAKE reference geometry built from the three bond lengths, and
  ``pq_input.in`` sets ``shake = mshake``. M-SHAKE solves the three coupled
  constraints of each water at once, so the default iteration counts are
  enough. Water without an H-H bond keeps two O-H constraints and a
  flexible angle.
- **1-4 scaling.** PQ has one Coulomb and one Lennard-Jones factor:
  :math:`1/\mathrm{SCEE}` and :math:`1/\mathrm{SCNB}` (AMBER),
  ``fudgeQQ`` and ``fudgeLJ`` (GROMACS), or the ratio of the CHARMM 1-4 to
  the normal well depth.

PQ computes Coulomb and Lennard-Jones with a shifted-force cutoff and no
long-range correction. Energies therefore differ from a PME or plain-cutoff
run of the same force field; the force-field terms themselves are the same.

Not supported
-------------

These terms stop the conversion with an error that names the term, because
PQ cannot evaluate them yet:

- CHARMM CMAP corrections, Urey-Bradley terms and harmonic impropers;
- non-uniform 1-4 scaling (for example GLYCAM) and atom-specific 1-4
  Lennard-Jones parameters;
- virtual sites and extra points (for example TIP4P water);
- AMOEBA terms, Ryckaert-Bellemans torsions and AMBER 10-12 hydrogen-bond
  terms.

The intra non-bonded file lists every pair inside a moltype, so it grows with
the square of the moltype size (about 45 000 pairs for a 300-atom protein).

Troubleshooting
---------------

``parm2pq needs ParmEd. Install it with 'pip install PQAnalysis[ff]' ...``
   ParmEd is missing. Run ``pip install "PQAnalysis[ff]"`` in the same
   environment.

``The force field contains 18 CHARMM CMAP correction maps; PQ cannot evaluate this term.``
   The force field uses a term from `Not supported`_; the message names it
   and counts it. AMBER ff19SB and the CHARMM protein force fields use CMAP.
   For proteins use AMBER ff14SB.

``The force field contains 1022 virtual sites or extra points (first: atom 308 'M'); ...``
   The water model has a massless site (TIP4P, TIP4P-Ew, OPC). Use a
   three-site model such as TIP3P or SPC/E.

``The force field has special 1-4 Lennard-Jones parameters ...`` or ``... uses non-uniform 1-4 scaling ...``
   PQ has one global 1-4 factor for Coulomb and one for Lennard-Jones.
   Force fields such as GLYCAM or CHARMM with pair-specific 1-4 values
   cannot be converted.

``AttributeError: 'NoneType' object has no attribute 'used'`` when ParmEd saves the prmtop
   The OpenMM system was built with ``constraints=HBonds`` or rigid water,
   so the constrained bonds have no force constant. Build it with
   ``constraints=None`` and ``rigidWater=False``; ``--shake-h`` adds the
   constraints in PQ. If the parameters go missing in another way,
   ``parm2pq`` itself stops with ``The bond 12-13 has no parameters.`` (or
   the same for an angle or dihedral).

``The structure has no coordinates; give a coordinate file.``
   A prmtop, top or psf has no positions. Pass them with ``-c``.

``parm2pq does not replace existing files without --overwrite: ...``
   The output directory already has setup files. Pick a new ``-d`` or add
   ``--overwrite``.

The box is too small, or the run has no water box
   Without a box in the coordinate file, ``parm2pq`` writes a vacuum box
   (extent plus two cutoffs). For a solvated system make sure the
   coordinate file has the box: an AMBER inpcrd from tleap or ParmEd does,
   a PDB without ``CRYST1`` does not. Each box length must be larger than
   twice ``rcoulomb``.

``intra_nonbonded.dat`` is very large
   The file lists every atom pair inside a moltype, so it grows with the
   square of the moltype size: about 45 000 pairs for a 300-atom protein,
   about 10\ :sup:`7` for a 5000-atom protein. Water and ions are separate
   small moltypes and add almost nothing. Large proteins work, but the file
   and the PQ setup time grow.

Validation
----------

Trp-cage (PDB 1L2Y, AMBER ff14SB) was converted from an AMBER prmtop and
compared with OpenMM, whose non-bonded forces were set to the PQ
shifted-force form with the same 12 Angstrom cutoff. The relative difference
of every energy term (bond, angle, dihedral, improper, Coulomb,
Lennard-Jones) was below :math:`10^{-5}` and the largest force difference
below :math:`10^{-6}` kcal mol\ :sup:`-1` Angstrom\ :sup:`-1`, in vacuum
and with 1249 TIP3P waters and one chloride ion. A copy of the vacuum
system with all torsion phases moved away from 0 and 180 degrees agreed as
well, which checks the sign of the phase.
