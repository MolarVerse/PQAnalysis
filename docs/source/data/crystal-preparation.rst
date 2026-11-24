Prepare a crystal from CIF
==========================

Use this workflow to turn an ordered periodic CIF into a PQ geometry restart.
It works for inorganic solids, molecular crystals, and ordered MOF frameworks.

.. code-block:: console

   pip install 'PQAnalysis[crystal]'
   pqanalysis crystal input.cif supercell.rst --repeat 2 2 1

This expands the CIF symmetry, repeats the cell along its first two lattice
vectors, and writes a centered PQ restart. A 20-atom expanded cell becomes
80 atoms. Without ``--repeat``, the command writes one expanded cell.
If the CIF contains multiple structure blocks, add ``--block NAME``, using
the block name without the ``data_`` prefix.

Before using the result
-----------------------

1. Open the restart in PQViewer and check the composition, atom count, and
   cell lengths and angles. Atoms on opposite faces may be periodic images.
2. Choose a compatible PQ mode and model. Replication preserves cell angles;
   it does not make a skewed cell orthogonal. Classical PQ MM-MD does not
   accept triclinic cells.
3. Check that the box suits the interaction cutoff. A successful export
   does not check cutoff compatibility or assign a topology, force field,
   or charges. To convert an already parameterized system, see
   :doc:`force-field-import`.

If a CIF is rejected
--------------------

For partial occupancies or equivalent source sites, inspect the original CIF
and resolve the sites explicitly, or use a curated ordered structure. Keep the
original and rerun the command on a copy. This tool does not remove solvent,
counterions, or framework atoms, choose a disorder model, or activate a MOF.
A structure without a valid periodic cell is rejected rather than padded with
vacuum.
