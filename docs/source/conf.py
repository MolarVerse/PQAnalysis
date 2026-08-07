"""Sphinx configuration for the PQAnalysis documentation."""

import sys

from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parent
DOCS_DIR = SOURCE_DIR.parent
PROJECT_ROOT = DOCS_DIR.parent

sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(SOURCE_DIR / "_plots"))

project = "PQAnalysis"
author = "the PQAnalysis authors"
copyright = "2023-2026, the PQAnalysis authors"

try:
    from PQAnalysis import __version__ as release
except Exception:  # pragma: no cover - package may be absent in a bare checkout
    release = ""
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.inheritance_diagram",
    "sphinx_sitemap",
    "matplotlib.sphinxext.plot_directive",
    "myst_parser",
    "sphinx_copybutton",
]

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

autoclass_content = "both"
autodoc_class_signature = "mixed"
autodoc_typehints_format = "short"
autodoc_member_order = "alphabetical"
maximum_signature_line_length = 50
add_module_names = False

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

plot_formats = [("svg", 96)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_include_source = False

templates_path = ["_templates"]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
exclude_patterns = []
highlight_language = "python"

html_theme = "furo"
html_title = "PQAnalysis"
html_logo = "logo/PQAnalysis.png"
html_favicon = "logo/PQAnalysis.png"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_baseurl = "https://molarverse.github.io/PQAnalysis/"

html_theme_options = {
    "sidebar_hide_name": False,
    "light_css_variables": {
        "color-brand-primary": "#1f718f",
        "color-brand-content": "#176c8c",
    },
    "dark_css_variables": {
        "color-brand-primary": "#65bddb",
        "color-brand-content": "#65bddb",
    },
    "source_repository": "https://github.com/MolarVerse/PQAnalysis/",
    "source_branch": "main",
    "source_directory": "docs/source/",
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/MolarVerse/PQAnalysis",
            "html": (
                '<svg stroke="currentColor" fill="currentColor" stroke-width="0" '
                'viewBox="0 0 16 16"><path fill-rule="evenodd" d="M8 0C3.58 0 0 '
                '3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01'
                '-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13'
                '-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 '
                '2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31'
                '-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 '
                '1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 '
                '1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 '
                '3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55'
                '.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>'
            ),
            "class": "",
        },
    ],
}


def show_inherited_mixins(app, what, name, obj, options, lines):
    """Show inherited mixins in the base classes of a class."""

    if what != "class" or not hasattr(obj, "__bases__"):
        return

    for base in obj.__bases__:
        if base.__name__.endswith("Mixin"):
            options["inherited-members"] = True


def run_apidoc(app):
    """Generate the complete package API reference."""
    import better_apidoc

    better_apidoc.APP = app
    better_apidoc.main([
        "better-apidoc",
        "-t",
        str(SOURCE_DIR / "_templates"),
        "--force",
        "--no-toc",
        "--separate",
        "-o",
        str(SOURCE_DIR / "code"),
        str(PROJECT_ROOT / "PQAnalysis"),
    ])


def setup(app):
    app.connect("autodoc-process-docstring", show_inherited_mixins)
    app.connect("builder-inited", run_apidoc)
