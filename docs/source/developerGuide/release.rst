Release Process
===============

PQAnalysis uses ``dev`` as the integration branch and releases from ``main``.
Version strings are derived from Git tags by ``setuptools_scm``; do not edit the
generated ``PQAnalysis/_version.py`` file.

Release boundary
----------------

The release workflow matches every pushed tag. A tag push can publish to PyPI
and TestPyPI, create a GitHub release, sign release artifacts and update
``CHANGELOG.md`` on ``main``.

.. warning::

   Treat ``git push origin vX.Y.Z`` as the publication action. Deleting a Git
   tag does not remove an uploaded Python distribution.

Pre-release checks
------------------

1. Open a release pull request from ``dev`` to ``main``.
2. Confirm that the release PR contains only reviewed integration changes.
3. Run both runtime-type-checking configurations with ``bash pytest.sh``.
4. Build the HTML documentation and link check with warnings as errors.
5. Confirm the Conventional Commits history produces meaningful release notes.
6. Verify the intended version is greater than every existing release tag.

The local verification commands are:

.. code-block:: console

   $ git fetch origin --tags
   $ bash pytest.sh
   $ python -m sphinx -E -W --keep-going \
       -b html docs/source docs/build/html
   $ python -m sphinx -E -W --keep-going \
       -b linkcheck docs/source docs/build/linkcheck

Tagging
-------

After the release PR is merged and the ``main`` checks pass, tag the verified
``main`` commit with the existing ``vMAJOR.MINOR.PATCH`` convention:

.. code-block:: console

   $ git switch main
   $ git pull --ff-only origin main
   $ git tag -a vX.Y.Z -m "PQAnalysis vX.Y.Z"
   $ git push origin vX.Y.Z

Use a major version for incompatible public API or file-contract changes, a
minor version for backward-compatible functionality and a patch version for
backward-compatible fixes.

Publication verification
------------------------

Do not consider the release complete until all of these are confirmed:

* the release workflow succeeded;
* the new version is available from PyPI;
* the GitHub release contains the signed distribution artifacts;
* ``CHANGELOG.md`` was updated on ``main``;
* the documentation workflow deployed the verified ``main`` build;
* a clean environment can install the published version and run
  ``pqanalysis --version``.

Published PyPI files are immutable. Correct a defective release with a new
patch release rather than moving or reusing its tag.
