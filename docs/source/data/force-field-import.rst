.. _forceFieldImport:

Force-Field Import (parm2pq)
============================

``parm2pq`` converts a parametrized system from AMBER, CHARMM or GROMACS into
the setup files of a PQ ``mm-md`` run. It reads the input with
`ParmEd <https://parmed.github.io/ParmEd/>`_, an optional dependency:

.. code-block:: console

   $ pip install PQAnalysis[ff]

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
