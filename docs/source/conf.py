# Configuration file for the Sphinx documentation builder.
#

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

import glob, re

# -- Project information -----------------------------------------------------

project = 'topiary'
copyright = '2021, Michael J. Harms'
author = 'Michael J. Harms'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
   'sphinx.ext.duration',
   'sphinx.ext.doctest',
   'sphinx.ext.autodoc',
   'sphinx.ext.autosummary',
   'sphinxcontrib.napoleon'
]

autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

#master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_logo = "_static/img/logo_150x150.png"

html_favicon = "_static/img/favicon.ico"

html_css_files = [
    'css/stylesheet.css',
]


# Files to exclude from final docs
exclude_patterns = ['_build', 'links.rst', 'topiary.rst','modules.rst']

# This allows us to put all of the links in a single .rst file
rst_epilog = ""
# Read link all targets from file
with open('links.rst') as f:
     rst_epilog += f.read()

# Make sure we actually run sphinx-apidoc on readthedocs
os.system("sphinx-apidoc -f ../../topiary -o .")

# This wackiness cleans up the sphinx-apidoc so it's much cleaner and easier
# to read
topiary_rst = glob.glob("../source/topiary.*.rst")

for rst in topiary_rst:

    header = []
    body = []
    footer = []

    wipe_out_next = False
    get_footer = False
    counter = 0
    with open(rst) as f:
        for line in f:

            if counter < 3:
                header.append(line)
                counter += 1
                continue

            if line.startswith("Submodules") or line.startswith("Subpackages"):
                wipe_out_next = True
                continue

            if wipe_out_next:
                wipe_out_next = False
                continue

            if line.startswith("Module contents"):
                get_footer = True

            if get_footer:
                footer.append(line)
                continue

            line = re.sub(" module","",line)
            body.append(line)

            #line = re.sub(" package","",line)
            #body.append(line)

    header[0] = re.sub(" package","",header[0])

    final_out = header[:]
    final_out.append("\n")
    final_out.extend(footer[2:])
    final_out.append("\n")
    final_out.extend(body)
    final_out.append("\n")

    g = open(rst,"w")
    g.writelines(final_out)
    g.close()
