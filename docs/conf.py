# Configuration file for the Sphinx documentation builder.
#
# For the full list of options see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = 'RABIES'
copyright = '2019, CoBrALab and Gabriel Desrosiers-Gregoire and Gabriel A. Devenyi and Mallar Chakravarty'
author = 'CoBrALab'

# The full version, including alpha/beta/rc tags
release = '0.6.1'

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx.ext.githubpages",
    "sphinx_rtd_dark_mode",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.todo",
    'sphinxcontrib.bibtex',
    'sphinxcontrib.programoutput',
]

# to get bibliography
bibtex_bibfiles = ['_static/refs.bib']

# Choose to generate TODO notices or not. Defaults to False
todo_include_todos = False

# Set MyST specific extensions
myst_enable_extensions = [
    "tasklist",
    "amsmath",
    "dollarmath",
    "colon_fence",  # ::: fences, so admonitions can nest and stay readable
    "deflist",
    "attrs_inline",
]

# enable equation rendering inline
myst_dmath_double_inline = True

# Generate anchors for h1-h3 so other pages can deep-link to a section
myst_heading_anchors = 3

# Make sure the target is unique
autosectionlabel_prefix_document = True

# Number and cross-reference figures, tables and code blocks
numfig = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------

# sphinx_rtd_dark_mode forces this theme, it is set here so the value is not
# silently disagreeing with what actually gets built.
html_theme = 'sphinx_rtd_theme'

html_title = f'RABIES {release}'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = ['custom.css']

# -- Options for sphinx_rtd_dark_mode ----------------------------------------
default_dark_mode = False

# -- Options for the RTD theme -----------------------------------------------
# The whole navigation structure is included in the sidebar of every page.
html_theme_options = {
    'collapse_navigation': False,
    'navigation_depth': 3,
    'sticky_navigation': True,
    'titles_only': False,
}

# -- Options for sphinx_copybutton -------------------------------------------
# Strip shell prompts and continuation markers when copying a command
copybutton_prompt_text = r'>>> |\.\.\. |\$ '
copybutton_prompt_is_regexp = True
copybutton_line_continuation_character = '\\'
