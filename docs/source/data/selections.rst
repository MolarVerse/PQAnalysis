Atom selections
===============

Analysis commands and Python objects select atoms with
:class:`PQAnalysis.topology.selection.Selection`. The same string language is
used in input-file keys such as ``reference_selection`` and
``target_selection``. Indices are **0-based**.

Always check that the selection contains the intended atoms. RDF normalization,
MSD statistics and VACF amplitudes all scale with that population.

String language
---------------

A selection string is parsed with a Lark grammar. Applied to
``examples/water/trajectory.xyz``, whose atoms are O (index 0), H (1) and
H (2), the common forms select:

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - String
     - Meaning
     - Atoms selected
   * - ``O``
     - Atom-type name
     - O
   * - ``H``
     - Atom-type name
     - both H
   * - ``0``
     - Single index
     - O
   * - ``1..2`` or ``1-2``
     - Inclusive index range
     - both H
   * - ``0-2``
     - Inclusive index range
     - all three
   * - ``0..2..2``
     - Range with step
     - O and the second H
   * - ``elem(8)`` or ``elem(O)``
     - Element number or symbol
     - O
   * - ``elem(H)``
     - Element symbol
     - both H
   * - ``atom(O, 8)``
     - Atom type and element
     - O
   * - ``*`` or ``all``
     - Every atom
     - all three
   * - ``*|H``
     - Set difference (all except H)
     - O
   * - ``O,H``
     - Union
     - all three
   * - ``* & elem(1)``
     - Intersection
     - both H

Operators
---------

Statements combine with:

``,``
   union

``&``
   intersection

``|``
   set difference (left minus right)

Use parentheses whenever you mix operators. The parser's implicit grouping
does not follow the order the class docstring describes (``O,H&H`` is read as
``(O,H)&H``; see `issue #185
<https://github.com/MolarVerse/PQAnalysis/issues/185>`_), so an
unparenthesized mix is not portable across versions.

Python
------

.. code-block:: python

   import numpy as np

   from PQAnalysis.io import read_trajectory
   from PQAnalysis.topology import Selection

   traj = read_trajectory("examples/water/trajectory.xyz")
   topology = traj[0].topology

   Selection("O").select(topology)          # array([0])
   Selection("H").select(topology)          # array([1, 2])
   Selection("*|H").select(topology)        # array([0])
   Selection(np.array([0, 2])).select(topology)  # array([0, 2])
   Selection(None).select(topology)         # all atoms

``Selection`` also accepts a single ``Atom`` or ``Element``. A plain Python
list of integers is rejected; pass a NumPy integer array.

``use_full_atom_info``
----------------------

By default ``select(..., use_full_atom_info=False)`` matches on element type
only. Set ``use_full_atom_info=True`` to distinguish atom-type names that share
an element (PQ residue atom types). Residue-aware RDF exclusions additionally
need a restart and a moldescriptor; see :doc:`../analyses/rdf`.

Analysis objects take the same strings as ``reference_species``,
``target_species`` or ``selection``:

.. code-block:: python

   from PQAnalysis.analysis import RDF
   from PQAnalysis.io import TrajectoryReader

   RDF(
       TrajectoryReader("examples/water/trajectory.xyz"),
       reference_species="O",
       target_species="H",
       delta_r=0.5,
       r_max=4.0,
   )
