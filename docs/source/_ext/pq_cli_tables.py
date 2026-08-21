"""
A small Sphinx extension that keeps the command-line reference in sync
with the code.

The ``pq-cli-table`` directive renders a list-table of commands. Every
command name in the table body is validated against the registry in
:py:mod:`PQAnalysis.cli.main` at build time, so a renamed or removed
command fails the documentation build instead of leaving a stale table
row. The ``pq-cli-covered`` directive marks commands that are documented
in prose instead of a table. After the build has read every page, the
extension checks that each registered command was documented exactly
once and raises a warning (an error under ``-W``) for every command
that is missing from or duplicated in the documentation.

Table rows use the form ``name -- purpose`` or
``name -- purpose -- primary input``; the purpose text is editorial,
while the command name and its cross reference come from the code.
"""

import ast
import contextlib
import functools
import importlib
import io

from pathlib import Path

import PQAnalysis.cli as cli_module

from docutils import nodes
from docutils.parsers.rst import directives
from docutils.statemachine import ViewList
from sphinx.util import logging
from sphinx.util.docutils import SphinxDirective

logger = logging.getLogger(__name__)

_ROW_SEPARATOR = " -- "


@functools.lru_cache(maxsize=1)
def _registered_commands():
    """
    Returns the registered command names from the PQAnalysis CLI.

    The dispatcher module is read with :py:mod:`ast` instead of being
    imported: during the documentation build the api-doc generator
    imports the package itself, so importing the dispatcher here can
    observe a partially initialized module whose class attributes do
    not exist yet. The class names of the dispatch table and their
    defining modules are therefore taken from the source, and only the
    individual command modules are imported to read their program
    names. Importing a command module prints the PQAnalysis header to
    stdout, which is swallowed so that it does not clutter the build
    output.

    Returns
    -------
    set[str]
        The names of all registered ``pqanalysis`` subcommands.

    Raises
    ------
    RuntimeError
        If no command could be read from the dispatcher module.
    """
    main_source = Path(cli_module.__file__).with_name("main.py")
    tree = ast.parse(main_source.read_text(encoding="utf-8"))

    modules_of_class = {}
    dispatched_classes = set()

    for node in ast.walk(tree):
        # 'from .rdf import RDFCLI' -> {'RDFCLI': 'rdf'}
        if isinstance(node, ast.ImportFrom) and node.module:
            for alias in node.names:
                modules_of_class[alias.asname or alias.name] = node.module

        # 'RDFCLI.program_name(): RDFCLI' inside the dispatch table
        if isinstance(node, ast.Attribute) and node.attr == "program_name":
            if isinstance(node.value, ast.Name):
                dispatched_classes.add(node.value.id)

    commands = set()

    with contextlib.redirect_stdout(io.StringIO()):
        for class_name in dispatched_classes:
            module_name = modules_of_class.get(class_name)

            if module_name is None:
                continue

            module = importlib.import_module(
                f"{cli_module.__name__}.{module_name}"
            )
            commands.add(getattr(module, class_name).program_name())

    if not commands:
        raise RuntimeError(
            "the pq-cli-table extension could not read any command from "
            f"{main_source}; the dispatch table format may have changed"
        )

    return commands


def _documented(env):
    """
    Returns the per-document map of documented command names.

    Parameters
    ----------
    env : sphinx.environment.BuildEnvironment
        The active build environment.

    Returns
    -------
    dict[str, list[str]]
        Command names documented per docname.
    """
    if not hasattr(env, "pq_cli_documented"):
        env.pq_cli_documented = {}

    return env.pq_cli_documented


class PQCliTable(SphinxDirective):
    """
    Renders a validated list-table of ``pqanalysis`` subcommands.
    """

    has_content = True
    required_arguments = 0
    optional_arguments = 0
    option_spec = {"title": directives.unchanged}

    def run(self):
        """
        Builds the list-table nodes from the directive content.

        Returns
        -------
        list[docutils.nodes.Node]
            The rendered table.
        """
        registry = _registered_commands()
        rows = []
        n_columns = 2

        for line in self.content:
            if not line.strip():
                continue

            parts = [part.strip() for part in line.split(_ROW_SEPARATOR)]

            if len(parts) not in (2, 3):
                raise self.error(
                    "pq-cli-table rows must be 'name -- purpose' or "
                    f"'name -- purpose -- input', got: {line!r}"
                )

            name = parts[0]

            if name not in registry:
                raise self.error(
                    f"pq-cli-table lists {name!r}, which is not a "
                    "registered pqanalysis command. Registered commands: "
                    f"{', '.join(sorted(registry))}"
                )

            _documented(self.env).setdefault(
                self.env.docname, []
            ).append(name)
            rows.append(parts)
            n_columns = max(n_columns, len(parts))

        title = self.options.get("title", "Commands")
        widths = "24 50 26" if n_columns == 3 else "28 72"
        headers = ["Command", "Purpose"]

        if n_columns == 3:
            headers.append("Primary input")

        text = ViewList()

        def emit(line):
            text.append(line, "pq-cli-table")

        emit(f".. list-table:: {title}")
        emit("   :class: pq-command-table")
        emit("   :header-rows: 1")
        emit(f"   :widths: {widths}")
        emit("")

        emit(f"   * - {headers[0]}")

        for header in headers[1:]:
            emit(f"     - {header}")

        for parts in rows:
            name = parts[0]
            padded = parts + [""] * (n_columns - len(parts))
            emit(f"   * - :ref:`{name} <cli.{name}>`")

            for cell in padded[1:]:
                emit(f"     - {cell}")

        node = nodes.section()
        node.document = self.state.document
        self.state.nested_parse(text, self.content_offset, node)

        return node.children


class PQCliCovered(SphinxDirective):
    """
    Marks commands as documented in prose instead of a table.

    The directive renders nothing; it only records its content so the
    completeness check accepts commands such as ``convert`` that are
    described in running text.
    """

    has_content = True

    def run(self):
        """
        Records the covered command names.

        Returns
        -------
        list
            An empty node list.
        """
        registry = _registered_commands()

        for line in self.content:
            name = line.strip()

            if not name:
                continue

            if name not in registry:
                raise self.error(
                    f"pq-cli-covered lists {name!r}, which is not a "
                    "registered pqanalysis command."
                )

            _documented(self.env).setdefault(
                self.env.docname, []
            ).append(name)

        return []


def _purge(app, env, docname):  # pylint: disable=unused-argument
    """
    Drops the recorded commands of a document that is re-read.
    """
    _documented(env).pop(docname, None)


def _merge(app, env, docnames, other):  # pylint: disable=unused-argument
    """
    Merges recorded commands from a parallel build worker.
    """
    for docname, names in _documented(other).items():
        _documented(env).setdefault(docname, []).extend(names)


def _check(app, env):
    """
    Verifies every registered command is documented exactly once.
    """
    documented = [
        name for names in _documented(env).values() for name in names
    ]

    if not documented:
        return

    registry = _registered_commands()
    missing = sorted(registry - set(documented))
    duplicated = sorted(
        {name for name in documented if documented.count(name) > 1}
    )

    if missing:
        logger.warning(
            "pqanalysis commands missing from the command-line "
            "reference: %s",
            ", ".join(missing),
        )

    if duplicated:
        logger.warning(
            "pqanalysis commands documented more than once in the "
            "command-line reference: %s",
            ", ".join(duplicated),
        )


def setup(app):
    """
    Registers the directives and consistency check.

    Parameters
    ----------
    app : sphinx.application.Sphinx
        The Sphinx application.

    Returns
    -------
    dict
        Extension metadata.
    """
    app.add_directive("pq-cli-table", PQCliTable)
    app.add_directive("pq-cli-covered", PQCliCovered)
    app.connect("env-purge-doc", _purge)
    app.connect("env-merge-info", _merge)
    app.connect("env-check-consistency", _check)

    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
