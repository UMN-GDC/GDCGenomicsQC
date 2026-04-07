# -- Project information -----------------------------------------------------

project = 'GDCgenomicsQC pipeline'
copyright = '2026, Genomic data commons'
author = 'Michael Anderson and Christian Coffman'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.mathjax',
    'sphinxcontrib.mermaid',
]

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

# Ensure MathJax uses a reliable CDN for rendering equations
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory.
html_static_path = ['_static']
