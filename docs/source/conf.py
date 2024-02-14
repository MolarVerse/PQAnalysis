# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
import os


sys.path.insert(0, os.path.abspath('../../'))

project = 'PQAnalysis'
copyright = '2023, Jakob Gamper, Josef M. Gallmetzer, Clarissa A. Seidler'
author = 'Jakob Gamper, Josef M. Gallmetzer, Clarissa A. Seidler'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx_sitemap',
    'sphinx.ext.inheritance_diagram',
    'myst_parser',
]

# Napoleon settings
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

autoclass_content = 'both'
autodoc_class_signature = 'mixed'
autodoc_typehints_format = 'short'
autodoc_member_order = 'alphabetical'
maximum_signature_line_length = 50
add_module_names = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ['.rst', '.md']

# The master toctree document.
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = []

highlight_language = 'python'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_style = 'css/custom.css'
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',  # Provided by Google in your dashboard
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,

    'logo_only': False,

    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'includehidden': True,
    'titles_only': True,
    'globaltoc_maxdepth': -1,
}

html_logo = 'logo/PQAnalysis.png'
# github_url = ''
# html_baseurl = ''

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


def show_inherited_mixins(app, what, name, obj, options, lines):
    """Show inherited mixins in the base classes of a class"""

    if what != 'class' or not hasattr(obj, '__bases__'):
        return

    for base in obj.__bases__:
        if base.__name__.endswith('Mixin'):
            options['inherited-members'] = True


def run_apidoc(app):
    """Generage API documentation"""
    import better_apidoc
    better_apidoc.APP = app
    better_apidoc.main([
        'better-apidoc',
        '-t',
        os.path.join('.', 'source', '_templates'),
        '--force',
        '--no-toc',
        '--separate',
        '-o',
        os.path.join('.', 'source', 'code'),
        os.path.join('..', 'PQAnalysis')
    ])


def setup(app):
    app.connect('autodoc-process-docstring', show_inherited_mixins)
    app.connect('builder-inited', run_apidoc)
