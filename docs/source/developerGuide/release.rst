Release Process
===============

PQAnalysis uses ``dev`` as the integration branch and releases from ``main``.
Version strings are derived from Git tags by ``setuptools_scm``; do not edit the
generated ``PQAnalysis/_version.py`` file.

Release boundary
----------------

The release workflow matches every pushed tag. A tag push can publish to PyPI,
create a GitHub release with the notes of that version and sign release
artifacts. It does not modify the repository: ``CHANGELOG.md`` and
``DEV-CHANGELOG.md`` are prepared in the release pull request.

.. warning::

   Treat ``git push origin vX.Y.Z`` as the publication action. Deleting a Git
   tag does not remove an uploaded Python distribution.

Pre-release checks
------------------

1. Merge ``main`` into ``dev`` when ``main`` has commits that ``dev`` lacks.
2. Run the *Create Release PR* workflow with the intended ``vX.Y.Z`` version.
   It creates ``release/vX.Y.Z`` from ``dev``, moves the changelog fragments
   into ``CHANGELOG.md`` and ``DEV-CHANGELOG.md`` and opens the pull request
   to ``main``.
3. Confirm that the release PR contains only reviewed integration changes and
   that the generated changelog entries are meaningful.
4. Run both runtime-type-checking configurations with ``bash pytest.sh``.
5. Build the HTML documentation and link check with warnings as errors.
6. Verify the intended version is greater than every existing release tag.

The local verification commands are:

.. code-block:: console

   $ git fetch origin --tags
   $ bash pytest.sh
   $ python -m sphinx -E -W --keep-going \
       -b html docs/source docs/build/html
   $ python -m sphinx -E -W --keep-going \
       -b linkcheck docs/source docs/build/linkcheck

Changelog fragments
-------------------

Every regular pull request adds a fragment under ``changes/user/`` or
``changes/developer/``; ``changes/README.md`` lists the categories. The
*Changelog Fragment Check* workflow enforces this on pull requests to ``dev``.
A maintainer can apply the ``skip-changelog`` label to a pull request that needs
no entry. The release pull request consumes the fragments, so ``CHANGELOG.md``
and ``DEV-CHANGELOG.md`` change only there.

The *Create Release PR* workflow pushes with ``GITHUB_TOKEN`` unless the
repository secret ``RELEASE_TOKEN`` is set. Pull requests and pushes made with
``GITHUB_TOKEN`` do not start other workflows, so set ``RELEASE_TOKEN`` to a
fine-grained token with contents and pull-request write access if the release
pull request should run its checks automatically.

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
* the GitHub release notes match the changelog entries of the release PR;
* the documentation workflow deployed the verified ``main`` build;
* a clean environment can install the published version and run
  ``pqanalysis --version``.

Published PyPI files are immutable. Correct a defective release with a new
patch release rather than moving or reusing its tag.
