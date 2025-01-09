# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from pathlib import Path

HERE = Path(__file__).parent
PROJECT_ROOT = HERE.parent
sys.path.insert(0, str(PROJECT_ROOT))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BurstLink'
copyright = '2025, LiyingZhou'
author = 'LiyingZhou'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'myst_nb',  
]
autodoc_mock_imports = [
    "sphinx",
    "sphinx_rtd_theme",
    "sphinx_mdinclude",
    "sphinx_autodoc_typehints",
    "sphinx_rtd_size",
    "nbsphinx",
    "sphinx_gallery",
    "furo",
    "myst_nb",
    "leidenalg",
    "matplotlib",
    "numpy",
    "pandas",
    "plotly",
    "POT",
    "pycirclize",
    "pydpc",
    "python_igraph",
    "scikit_learn",
    "scipy",
    "seaborn",
    "tqdm",
    "gseapy",
]
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
