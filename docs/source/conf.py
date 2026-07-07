import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

project = "MOT-MIC"
author = "yuhan1li"
copyright = "2026, yuhan1li"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "myst_nb",
]

autosummary_generate = True
templates_path = ["_templates"]
exclude_patterns = []

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "_static/img/motmic_head_image.png"

nb_execution_mode = "off"

