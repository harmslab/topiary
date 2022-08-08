
.. include:: links.rst

.. role:: emph

.. _drawing-trees-doc:

=============
Drawing Trees
=============

Topiary allows users to generate their own plots with customized formatting
using the topiary python :ref:`API<api-doc>`. Users can control all aspects of
phylogenetic tree drawing from a the :code:`topiary.draw.tree` function. The
following python generates the following output

.. code-block:: python

  import topiary
  topiary.draw.tree("06_ancestors_with_branch_supports")

.. image:: _static/img/final-tree.svg
  :align: center
  :alt: Topiary tree drawing

.. note::

  Jupyter examples coming soon.


Description of arguments
========================

---------------------
General plot features
---------------------

+ :code:`output_file`. Output file. Type is determined by file extension.
  Allowed extensions are pdf, svg, and png. If this is not specified and you are
  running in a jupyter notebook, the tree will be written out to the screen only.
+ :code:`font_size`. Controls overall font size for the document. This is the
  size assigned to the tip labels. Node labels are given a 0.75 this size
  unless a specific label size is given using :code:`label_text_style`.
+ :code:`stroke_width`. Stroke width, in pixels, for the tree.
+ :code:`vertical_pixels_per_tip`. How many vertical pixels to assign each
  tip on the tree. This controls the height of the plot.
+ :code:`min_height`. Make sure the plot is at least this many pixels high.

--------------------------
Internal nodes (ancestors)
--------------------------

Labels
------

+ :emph:`What labels to draw?`

  + :code:`bs_label`. Draw bootstrap labels (:code:`True` or :code:`False`)
  + :code:`pp_label`. Draw posterior probability labels (:code:`True` or
    :code:`False`)
  + :code:`event_label`. Draw reconciliation event labels (:code:`True` or
    :code:`False`)
  + :code:`anc_label`. Draw ancestral name labels (:code:`True` or
    :code:`False`)

+ :emph:`How to format labels?`

  + :code:`label_position`. Where to position the labels relative to the node.
    ("top-left", "top", "top-right", "right", "bottom-right", "bottom",
    "bottom-left", or "left".)
  + :code:`label_position_offset`. How far, in pixels, to shift the position
    of the label in the direction indicated by label position.
  + :code:`label_color`. Label color. For details on how to specify color,
    see the :ref:`colors section<colors section>` below.
  + :code:`label_text_style`. This gives detailed control over the node label
    styling. This is a dictionary that contains css properties for the text.
    For details see the :ref:`CSS section<CSS Section>` below.

Colors
------

+ :code:`bs_color`. Branch support color map. Should be a dictionary with two
  elements that keys bootstrap values to color. For example,
  :code:`bs_color={50:"white",100:"black"}` would create a color gradient
  between white and black for bootstrap values between 50 and 100. Values
  above and below the input bootstrap values are set to the color of the
  relevant input value. In this case, a bootstrap value of 20 would be
  colored white. For details on how to specify color, see the
  :ref:`colors section<colors section>` below.
+ :code:`pp_color`. Posterior probability color map. See description of
  :code:`bs_color` for format.
+ :code:`event_color`. Color nodes according to the reconciliation event. This
  should be a dictionary keying events to colors. Allowable events are
  speciation (:code:`S`), duplication (:code:`D`), loss (:code:`L`), and
  lateral transfer (:code:`T`). For example, :code:`event_color={"D":"pink"}`
  would color all nodes that are duplications pink. Events that have no key
  in the dictionary will not be colored. For details on how to specify color,
  see the :ref:`colors section<colors section>` below.
+ :code:`node_color`. Set single color to all nodes. This overrides all other
  color arguments. For details on how to specify color, see the
  :ref:`colors section<colors section>` below.

Size
----

+ :code:`node_size`. Set size of nodes in pixels

---------------------------
Tree tips (modern proteins)
---------------------------

:emph:`How to construct names?`

+ :code:`tip_columns`. Which columns from the dataframe to use to build tip
  the tip labels.
+ :code:`tip_name_separator`. What to place between the tip column values.
  This can be any character (i.e. :code:`"|"`, :code:`","`, etc.).
+ :code:`disambiguate_tip_names`. Whether or not to append uid to duplicate
  tip names to make them unique (:code:`True` or :code:`False`).

The defaults are :code:`tip_columns=["species","recip_paralog"]` and
:code:`tip_name_separator="|"`, which gives tip names like "Homo sapiens|LY96".

:emph:`How to format tip labels?`

+ :code:`tip_text_style`. This gives detailed control over the tip label
  styling. This is a dictionary that contains css properties for the text.
  For details see the :ref:`CSS section<CSS section>` below.

Even more customization
=======================

The :code:`topiary.draw.tree` object is a wrapper that constructs a
:code:`topiary.draw.PrettyTree` object. Advanced users can create custom
instances of this object to further customize tree outputs. Please see the
:ref:`API<api-doc>` and :code:`topiary.draw.PrettyTree` docstring for details.

Topiary draws trees using the `toytree <toytree-link_>`_ library, which is
built on the `toyplot <toyplot-link_>`_ library. The :code:`topiary.draw.tree`
function returns a toyplot `canvas object <toyplot-canvas-link_>`_. Advanced
users can thus further manipulate the plot by editing that canvas as described
in the `toyplot <toyplot-link_>`_ and `toytree <toytree-link_>`_ documentation.


General considerations
======================

Because topiary trees are ultimately built using `toyplot <toyplot-link_>`_,
any customizations of :ref:`color<colors section>` or styling via
:ref:`css <CSS section>` accepted by toyplot are available to topiary users. For
full details, see the `toyplot documentation <toyplot-link_>`_. We summarize
the most commonly used options here.

.. _colors section:

------
Colors
------

Colors can be specified in three different ways:

+ RGBA tuples. Examples: :code:`(255,255,255,255)` would be white;
  :code:`(255,0,0,255)` would be red; :code:`(0,255,0,125)` would be half-opaque
  green.
+ Named colors. To get a full list of the colors available, you can run the
  following code in a jupyter notebook or python session.

  .. code-block:: python

    import toyplot
    print(toyplot.color.css.names)

+ Hexadecimal strings. Examples include :code:`"#407E98"` (a slate blue),
  :code:`"#FFFFFF"` (white), :code:`"#000000"` (black).

For a complete discussion of the color options, see the
`toyplot documentation <toyplot-colors-link_>`_.

.. _CSS section:

---
CSS
---

`CSS <css-wikipedia_>`_ (Cascading Style Sheets) is the language used to format
html in web pages. The `toyplot library <toyplot-link_>`_, that topiary uses
under the hood can a fair number of css properties and apply them to text. For
an up-to-date list of available properties, see the
`toyplot documentation <toyplot-css-link_>`_.
As of toyplot 1.02, the following css properties can be set:

+ "alignment-baseline"
+ "baseline-shift"
+ "fill"
+ "fill-opacity"
+ "font-family"
+ "font-size"
+ "font-weight"
+ "line-height"
+ "opacity"
+ "stroke"
+ "stroke-opacity"
+ "stroke-width"
+ "text-anchor"
+ "text-decoration-line"
+ "text-shadow"

For information on each of these css properties, there are css references online
(i.e. `w3schools <css-w3s_>`_).
