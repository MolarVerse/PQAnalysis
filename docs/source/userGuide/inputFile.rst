.. _inputFile:

Analysis Input Files
====================

RDF, MSD, VACF and vibrational analyses use a compact key-value format parsed
with `Lark <https://lark-parser.readthedocs.io/>`_. Each analysis documents its
required and optional keys in the generated command and input-reader reference.

Inline statements
-----------------

An inline statement assigns one value to one key:

.. code-block:: text

   key = value

Several assignments may share a line when separated by commas:

.. code-block:: text

   window = 1000, gap = 10, time_step = 0.001

Use separate lines for scientific input files unless a compact generated file
is required; one assignment per line is easier to review and diff.

Multiline lists
---------------

A key followed by values on subsequent lines creates a list. Terminate the
list with ``END``:

.. code-block:: text

   traj_files
   trajectory-001.xyz
   trajectory-002.xyz
   trajectory-003.xyz
   END

Multiline syntax always produces a list-like value. Use inline syntax for keys
that accept only a scalar.

Comments and case
-----------------

Keys are case-insensitive. The closing ``END`` token must be uppercase. Values
remain case-sensitive because they may contain filenames or selection
expressions. Leading and trailing whitespace is ignored. ``#`` starts a
comment, including at the end of a statement:

.. code-block:: text

   target_selection = O  # oxygen atoms

Value conversion
----------------

The parser converts strings to the type required by each documented key.
Common forms include:

.. list-table:: Input value forms
   :header-rows: 1
   :widths: 30 30 40

   * - Input
     - Parsed form
     - Notes
   * - ``True`` or ``False``
     - Boolean
     - Case-insensitive
   * - ``1``
     - Integer
     - May also satisfy a real-valued key
   * - ``1.0``
     - Floating-point number
     - Scientific notation is accepted where numeric keys permit it
   * - ``[1, 2, 3]``
     - List
     - Elements must have compatible types
   * - ``1..4`` or ``1-4``
     - ``range(1, 4)``
     - The stop value follows Python's exclusive convention
   * - ``1..3..10`` or ``1-3-10``
     - ``range(1, 10, 3)``
     - The middle value is the step
   * - ``frame-*.xyz``
     - Matching file list
     - Expanded with Python glob semantics

Filenames
---------

Relative filenames are interpreted from the command's working directory. The
ordinary filename grammar accepts letters, digits, ``_``, ``-`` and ``.``;
``*`` provides glob matching. For a portable analysis directory, keep the
input file and its referenced data together and run the command from that
directory.

Complete example
----------------

.. code-block:: text

   # oxygen-hydrogen radial distribution
   traj_files
   run-001.xyz
   run-002.xyz
   END

   reference_selection = O
   target_selection = H
   delta_r = 0.05
   r_max = 8.0
   out_file = rdf.dat

Run it with:

.. code-block:: console

   $ pqanalysis rdf rdf.in

See :doc:`../analyses/index` for analysis-specific examples and
:doc:`analysisOutputFiles` for output formats and scientific schemas.
