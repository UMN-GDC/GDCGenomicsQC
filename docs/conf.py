# -- Project information -----------------------------------------------------

project = 'GDCgenomicsQC pipeline'
copyright = '2026, Genomic data commons'
author = 'Michael Anderson and Christian Coffman'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.mathjax',
    'sphinxcontrib.mermaid',
    'myst_parser',
    'sphinx_design',
]

# The master toctree document.
root_doc = master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

# Ensure MathJax uses a reliable CDN for rendering equations
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"

# MyST Markdown configuration
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_image",
]
myst_heading_anchors = 3

# html_static_path = ['_static']  # disabled: _static directory does not exist
