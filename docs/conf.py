# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'RABIES Documentation'
copyright = '2019, CoBrALab and Gabriel Desrosiers-Gregoire and Gabriel A. Devenyi and Mallar Chakravarty'
author = 'CoBrALab'

# The full version, including alpha/beta/rc tags
release = '0.5.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx.ext.githubpages",
    "sphinx_rtd_dark_mode",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.todo",
    'sphinxcontrib.bibtex',
    'sphinxcontrib.programoutput',
]

# to get bibliography
bibtex_bibfiles = ['_static/refs.bib']

# Choose to generate TOOD notices or not. Defaults to False
todo_include_todos = False

# Set MyST specific extensions
myst_enable_extensions = [
    "tasklist",
    "amsmath",
    "dollarmath",
]

# enable equation rendering inline
myst_dmath_double_inline = True

# Make sure the target is unique
autosectionlabel_prefix_document = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'groundwork'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Options for sphinx_rtd_dark_mode -------
default_dark_mode = False

# Set some RTD theme config.  This includes the entire navigation structure
# into the sidebar of all pages.  However, expanding the sections isn't
# provided yet on the RTD theme (see
# https://github.com/readthedocs/sphinx_rtd_theme/issues/455).
html_theme_options = {
    'collapse_navigation': False,
    'navigation_depth': 3,
}

# Set some RTD theme config.  This includes the entire navigation structure
# into the sidebar of all pages.  However, expanding the sections isn't
# provided yet on the RTD theme (see
# https://github.com/readthedocs/sphinx_rtd_theme/issues/455).
html_theme_options = {
    'collapse_navigation': False,
    'navigation_depth': 2,
}
