# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent)]

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BurstLink'
copyright = '2024, LiyingZhou'
author = 'LiyingZhou'
release = '0.0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'myst_nb',  
]
autosummary_generate = True
nb_execution_mode = "off"  
source_suffix = {
    '.rst': 'restructuredtext',
    '.ipynb': 'myst-nb',
}
templates_path = ['_templates']
exclude_patterns = [
    'build',
    '_build',
    '_build/jupyter_execute/**',
    'Thumbs.db',
    '.DS_Store',
    'api/burstlink.rst', 
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = "furo"  
html_static_path = ["_static"]
html_css_files = [
    "css/override.css", 
]

html_logo = "_static/image/logo.png"

html_theme_options = {
    "sidebar_hide_name": True,
    "light_css_variables": {
        "color-brand-primary": "#357473",
        "color-brand-content": "#357473",
    },
}
