# Configuration file for the Sphinx documentation builder.
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from pathlib import Path

# 设置模块路径
HERE = Path(__file__).parent
PROJECT_ROOT = HERE.parent
MODULE_PATH = PROJECT_ROOT / "burstlink"  # 假设模块在项目根目录的 burstlink 文件夹中
if MODULE_PATH.exists():
    sys.path.insert(0, str(MODULE_PATH.parent))  # 将模块的父目录添加到 sys.path

# 验证模块是否能正确导入
try:
    import burstlink
    import burstlink.plotting
    import burstlink.preprocessing
    import burstlink.tools
    print("Modules imported successfully!")
except ImportError as e:
    print(f"Error importing module: {e}")

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BurstLink'
copyright = '2024, LiyingZhou'
author = 'LiyingZhou'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'myst_nb',  # 支持 Jupyter Notebook 文件
]
nb_execution_mode = "off"  # 禁用 Notebook 执行
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

# -- Autosummary configuration -----------------------------------------------
autosummary_generate = True  # 自动生成 .rst 文件
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': False,
    'show-inheritance': True,
}

# 添加的路径配置
sys.path.append(str(PROJECT_ROOT))
