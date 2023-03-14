"""
Class for drawing formatted phylogenetic trees using toytree.
"""

import topiary
from topiary.draw.core import construct_colormap
from topiary.draw.core import construct_sizemap
from topiary.draw.core import get_round_to
from topiary.draw.core import ete3_to_toytree
from topiary.draw.core import parse_position_string
from topiary.draw.core import color_to_css
import topiary._private.check as check

import toytree
import toyplot
from toyplot import svg
from toyplot import pdf

import numpy as np

import re
import copy
import xml


class PrettyTree:
    """
    Class for drawing formatted phylogenetic trees using toytree.

    Parameters
    ----------
    T : ete3.Tree or dp.Tree or toytree.tree or newick
        tree to draw. We *strongly* recommend this tree have uid as its
        tip labels to avoid mangled trees. If you want to assign prettier
        names to the tips, pass in name_dict.
    name_dict : dict, optional
        dictionary mapping strings in node.name to more useful names.
    font_size : float, default=15
        font size (pixels)
    stroke_width : float, default=2
        width of lines drawing tree (pixels)
    vertical_pixels_per_tip : float, default=20
        number of pixels to assign to each tip when calculating figure
        height
    aspect_ratio_parameters : tuple, default=(0.02,0.47)
        these parameters used to relate figure aspect ratio to number of
        taxa by the following equation:
        :code:`aspect_ratio_parameters[0]*num_taxa + aspect_ratio_parameters[1]`
    min_height : float, default=300
        minimum height for figure (pixels)
    height : float, optional
        height in pixels. If specified, overrides vertical_pixels_per_tip,
        aspect_ratio_parameters, and min_height.
    width : float, optional
        width in pixels. if specified, overrides aspect_ratio_parameters
    padding : float, optional
        padding in pixels.
    edge_style : dict, optional
        dictionary of css key/value pairs defining how the tree edge style is
        drawn
    draw_debug_points : bool, default=False
        draw points on the plot that validate relative x-scale of two plot
        areas.
    **kwargs : dict
        pass any other keyword arguments directly to toytree.tree.draw.
    """

    def __init__(self,
                 T,
                 name_dict=None,
                 font_size=15,
                 stroke_width=2,
                 vertical_pixels_per_tip=20,
                 aspect_ratio_parameters=(0.02,0.47),
                 min_height=300,
                 height=None,
                 width=None,
                 padding=None,
                 edge_style=None,
                 tip_labels_style=None,
                 draw_debug_points=False,
                 **kwargs):
        """
        Initialize PrettyTree class.
        """

        # ----------------------------------------------------------------------
        # Load tree and give pretty tip names

        # Convert tree into a toytree tree from ete3,
        if issubclass(type(T),toytree.tree):
            self._tT = T.copy()
        else:
            # Convert to a toytree
            self._tT = ete3_to_toytree(topiary.io.read_tree(T))

        # Rename tree.name entries according to name_dict
        if name_dict is not None:
            if not issubclass(type(name_dict),dict):
                err = "name_dict shouold be a dictionary mapping uid to pretty\n"
                err += "tip names.\n"
                raise ValueError(err)

            for idx in self._tT.idx_dict:
                name = self._tT.idx_dict[idx].name
                try:
                    self._tT.idx_dict[idx].name = name_dict[name]
                except KeyError:
                    pass

        # ----------------------------------------------------------------------
        # Artwork parameters

        self._font_size = check.check_float(font_size,
                                            "font_size",
                                            minimum_allowed=0)
        self._stroke_width = check.check_float(stroke_width,
                                               "stroke_width",
                                               minimum_allowed=0)
        self._vertical_pixels_per_tip = check.check_float(vertical_pixels_per_tip,
                                                            "vertical_pixels_per_tip",
                                                            minimum_allowed=0)

        aspect_ratio_parameters = check.check_iter(aspect_ratio_parameters,
                                                   "aspect_ratio_parameters",
                                                   minimum_allowed=2,
                                                   maximum_allowed=2)
        aspect_ratio_parameters = list(aspect_ratio_parameters)
        aspect_ratio_parameters[0] = check.check_float(aspect_ratio_parameters[0],
                                                       "aspect_ratio_parameters[0]")
        aspect_ratio_parameters[1] = check.check_float(aspect_ratio_parameters[1],
                                                       "aspect_ratio_parameters[1]")

        min_height = check.check_float(min_height,
                                       "min_height",
                                       minimum_allowed=0)

        # ----------------------------------------------------------------------
        # If height not passed directly, figure out

        if height is None:
            num_nodes = self._tT.ntips
            self._height = num_nodes*self._vertical_pixels_per_tip
            if self._height < min_height:
                self._height = min_height

        else:
            self._height = check.check_float(height,
                                             "height",
                                             minimum_allowed=0,
                                             minimum_inclusive=False)

        # ----------------------------------------------------------------------
        # If width not passed directly, figure out

        if width is None:
            ratio = aspect_ratio_parameters[0]*num_nodes + aspect_ratio_parameters[1]
            self._width = self._height/ratio
        else:
            self._width = check.check_float(width,
                                            "width",
                                            minimum_allowed=0,
                                            minimum_inclusive=False)

        # ----------------------------------------------------------------------
        # If padding not passsed directly, figure out

        if padding is None:
            self._padding = self._vertical_pixels_per_tip
        else:
            self._padding = check.check_float(padding,
                                              "padding",
                                              minimum_allowed=0)

        # ----------------------------------------------------------------------
        # set stroke and stroke-width in edge_style

        if edge_style is None:
            edge_style = {}

        if not issubclass(type(edge_style),dict):
            err = "edge_style should be a dictionary of css key/value pairs\n"
            raise ValueError(err)

        if "stroke_width" not in edge_style:
            edge_style["stroke-width"] = self._stroke_width

        if "stroke" not in edge_style:
            edge_style["stroke"] = "black"

        kwargs["edge_style"] = edge_style

        # ----------------------------------------------------------------------
        # set tip_labels_style

        if tip_labels_style is None:
            tip_labels_style = {}

        if not issubclass(type(tip_labels_style),dict):
            err = "tip_labels_style should be a dictionary of css key/value pairs\n"
            raise ValueError(err)

        if "font-size" not in tip_labels_style:
            tip_labels_style["font-size"] = f"{self._font_size}px"

        kwargs["tip_labels_style"] = tip_labels_style

        # ----------------------------------------------------------------------
        # Create canvas and axes on which to draw tree

        # figure out legend dimensions
        legend_height = min_height*0.167

        # Create canvas
        total_height = self._height + 3*self._padding + legend_height
        total_width = self._width + 2*self._padding
        self._canvas = toyplot.Canvas(height=f"{total_height}px",
                                      width=f"{total_width}px")

        # Create tree axis
        x1 = self._padding
        x2 = self._width - self._padding
        y1 = self._padding
        y2 = self._height - self._padding/2

        tree_del_x_range = x2 - x1
        tree_del_y_range = y2 - y1

        # Create toyplot axes for tree itself
        self._tree_ax = self._canvas.cartesian(bounds=(x1,x2,y1,y2))
        self._tree_ax.x.ticks.show = False
        self._tree_ax.y.ticks.show = False
        self._tree_ax.show = False

        # ----------------------------------------------------------------------
        # Draw tree skeleton and make preliminary render so we can calculate
        # conversions between domain (tree) coordinates and range (pixel)
        # coordinates

        self._draw_kwargs = kwargs
        _, self._tree_ax, self._tree_mark = self._tT.draw(axes=self._tree_ax,
                                                           **self._draw_kwargs)
        _ = toyplot.html.render(self.canvas)

        # Get domain (tree) to range (px) interconversion factors
        x_domain = np.array(self._tree_mark.domain("x"))
        x_range = self._tree_ax.project("x",x_domain)
        self._x_tree_to_px = np.abs((x_range[1] - x_range[0])/(x_domain[1] - x_domain[0]))
        self._x_px_to_tree = 1/self._x_tree_to_px

        y_domain = np.array(self._tree_mark.domain("y"))
        y_range = self._tree_ax.project("y",y_domain)
        self._y_tree_to_px = np.abs((y_range[1] - y_range[0])/(y_domain[1] - y_domain[0]))
        self._y_px_to_tree = 1/self._y_tree_to_px

        # By number of pixels we gave the tree (tree_del_x_range) and
        # conversion factor (_x_range_to_domain) we can calculate the dimensions
        # of the plot in tree space.
        self._x_total_domain = tree_del_x_range*self._x_px_to_tree
        self._x_min = x_domain[0]
        self._x_max = self._x_min + self._x_total_domain

        # Same for y
        self._y_total_domain = tree_del_y_range*self._y_px_to_tree
        self._y_min = y_domain[0]
        self._y_max = self._y_min + self._y_total_domain

        # Natural step in y in tree coordinates
        self._y_step = self._font_size*1.333*self._y_px_to_tree
        self._pixel_aspect = self._x_tree_to_px/self._y_tree_to_px

        # ----------------------------------------------------------------------
        # Create axes on which to draw legend

        # Since we set legened to be 0 to 10, 10 is top
        self._legend_y = 10
        self._legend_width = (x_domain[0]+self._x_total_domain)/2

        x1 = self._padding
        x2 = self._width - self._padding
        y1 = self._height + self._padding/2
        y2 = y1 + legend_height
        self._legend_ax = self._canvas.cartesian(bounds=(x1,x2,y1,y2),
                                                 xmin=x_domain[0],
                                                 xmax=x_domain[0]+self._x_total_domain,
                                                 ymin=0,
                                                 ymax=self._legend_y)

        self._legend_ax.x.ticks.show = False
        self._legend_ax.y.ticks.show = False
        self._legend_ax.show = False

        if draw_debug_points:
            self._tree_ax.scatterplot((x_domain[0],x_domain[1]),(0,0))
            self._legend_ax.scatterplot((x_domain[0],x_domain[1]),(10,10))


        # Dictionary to hold named nodes added
        self._plotted_nodes = {}


    def _get_node_values(self,
                         get_ancestors,
                         get_leaves,
                         get_root=True,
                         property_labels=None):
        """
        Get the x,y coordinates and properties of nodes.

        Parameters
        ----------
        get_ancestors : bool
            whether or not to get ancestors
        get_leaves : bool
            whether or not to get leaves
        get_root : bool, default=True
            whether or not to get the root node
        property_labels : str or list-like
            list of property labels to grab, in order. can take a single string
            as equivalent to a length-1 list.

        Returns
        -------
        x : numpy.ndarray
            x coordinates of nodes
        y : numpy.ndarray
            y coordinates of nodes
        all_props : dict or None
            dictionary keying property_labels to an np.ndarray of that property
            value. if property_labels is None, return None.
        """

        # Create mask with nodes to get
        idx = np.array(list(self._tT.idx_dict.keys()),dtype=int)
        leaf_mask = []
        root_mask = []
        for k in self._tT.idx_dict:
            leaf_mask.append(self._tT.idx_dict[k].is_leaf())
            root_mask.append(self._tT.idx_dict[k].is_root())

        leaf_mask = np.array(leaf_mask,dtype=bool)
        anc_mask = np.logical_not(leaf_mask)

        if get_ancestors:
            if get_leaves:
                mask = np.logical_or(leaf_mask,anc_mask)
            else:
                mask = anc_mask
        else:
            if get_leaves:
                mask = leaf_mask
            else:
                err = "\nplot_ancestors and/or plot_leaves must be True\n\n"
                raise ValueError(err)

        # If we are not getting root node
        if not get_root:
            root_mask = np.logical_not(np.array(root_mask,dtype=bool))
            mask = np.logical_and(root_mask,mask)

        # Get x/y coord for nodes
        node_coord = self._tT.get_node_coordinates()

        x = node_coord[idx[mask],0]
        y = node_coord[idx[mask],1]

        # Get property(s) of interest
        all_props = {}
        if property_labels is not None:

            # If passed in as a single string, put in a list to allow iteration
            # over single label
            if issubclass(type(property_labels),str):
                property_labels = [property_labels]

            # Iterate over all labels
            for p_label in property_labels:

                prop = []
                for i in idx[mask]:

                    node = self._tT.idx_dict[i]
                    if p_label not in node.features:
                        prop.append(None)
                        continue

                    try:
                        prop.append(self._tT.idx_dict[i].__dict__[f"_{p_label}"])
                    except KeyError:
                        prop.append(self._tT.idx_dict[i].__dict__[f"{p_label}"])

                all_props[p_label] = np.array(prop)

        else:
            all_props = None

        return x, y, all_props


    def draw_nodes(self,
                   property_label=None,
                   color=None,
                   size=None,
                   prop_span=None,
                   plot_ancestors=True,
                   plot_leaves=False,
                   plot_root=True,
                   scatter_style=None,
                   palette=None):
        """
        Draw nodes colored and/or sized by a property loaded into the tree.

        Parameters
        ----------
        property_label : str, optional
            color/size nodes by property_label, where property_label is a
            feature of nodes in the toytree.tree instance.
        color : str or tuple or dict
            set node color. If a single value, color all nodes that color. If
            list-like and length 2, treat as colors for minimum and maximum of a
            color gradient.  If dict, map property keys to color values. Colors
            can be RGBA tuples, named colors, or hexadecimal strings. See the
            toyplot documentation for allowed values.
        size : float or tuple or dict, optional
            set node size. If a single value, make all nodes that size. If
            list-like and length 2, treat as sizes for minimum and maximum of a
            size gradient. If dict, map property keys to size values. Sizes must
            be float >= 0.
        prop_span : tuple, optional
            set min/max values for property color/size calculation. First element
            is min, second is max. Only applicable if color and size are set
            to be gradients.
        plot_ancestors : bool, default=True
            draw ancestral nodes
        plot_leaves : bool, default=False
            draw leaf nodes
        plot_root : bool, default=True
            draw the root node
        scatter_style : dict, optional
            dictionary specifying how to draw the nodes. valid keys are "shape",
            "mstyle","angle","label", and "lstyle". see toyplot Markers
            documentation for description of values.
        palette : toyplot.color.Palette, optional
            custom toyplot palette. if specified, use this palette rather than
            color_span to define color scheme.
        """

        # Check property_label arg
        if property_label is not None:
            if not issubclass(type(property_label),str):
                err = "\nproperty_subclass should be a string pointing to a\n"
                err += "feature in the toytree.tree instance.\n\n"
                raise ValueError(err)

        # Check plot args
        plot_ancestors = check.check_bool(plot_ancestors,"plot_ancestors")
        plot_leaves = check.check_bool(plot_leaves,"plot_leaves")
        plot_root = check.check_bool(plot_root,"plot_root")

        # Get x, y, and property values
        x, y, prop_dict = self._get_node_values(plot_ancestors,
                                                plot_leaves,
                                                plot_root,
                                                property_label)

        # Get property from prop_dict
        if prop_dict is None:
            prop = np.ones(len(x),dtype=float)
        else:
            prop = prop_dict[property_label]

        # If no nodes, return self
        if len(prop) == 0:
            return self

        # Do not plot anything with "None" entry
        good_mask = np.array([p is not None for p in prop],dtype=bool)
        x = x[good_mask]
        y = y[good_mask]
        prop = prop[good_mask]

        # Don't plot if there are no nodes to plot
        if len(x) == 0:
            return self

        # Check the property span arg
        if prop_span is not None:

            prop_span = check.check_iter(prop_span,
                                         "prop_span",
                                         minimum_allowed=2,
                                         maximum_allowed=2)
            prop_span = list(prop_span)
            prop_span[0] = check.check_float(prop_span[0],
                                             "prop_span[0]")
            prop_span[1] = check.check_float(prop_span[1],
                                             "prop_span[1]")

            try:
                prop = np.array(prop,dtype=float)
            except ValueError:
                err = "\nproperty must be coercable to a float if prop_span is\n"
                err += "specified\n"
                raise ValueError(err)

            # Cap prop below min and above max
            prop[prop < prop_span[0]] = prop_span[0]
            prop[prop > prop_span[1]] = prop_span[1]

        # Construct colormap
        if color is None:
            color = "gray"

        # If color is a dictionary, drop any nodes that do not have that key
        if issubclass(type(color),dict):
            keys = list(color.keys())
            good_mask = np.array([p in keys for p in prop],dtype=bool)
            x = x[good_mask]
            y = y[good_mask]
            prop = prop[good_mask]

        cm, cm_span = construct_colormap(color,prop,prop_span,palette)

        # Construct size map
        if size is None:
            size = self.default_size

        sm, sm_span = construct_sizemap(size,prop,prop_span)

        # Construct size and color lists for each node. If there is no color or
        # size for this value, drop them from the plot
        sizes = []
        colors = []
        for p in prop:
            try:
                sizes.append(sm(p))
            except (KeyError,ValueError):
                sizes.append(None)

            try:
                colors.append(cm(p))
            except (KeyError,ValueError):
                colors.append(None)

        keep_mask = np.logical_and(np.array([s is not None for s in sizes],dtype=bool),
                                   np.array([c is not None for c in colors],dtype=bool))

        x = x[keep_mask]
        y = y[keep_mask]
        sizes = np.array(sizes)[keep_mask]
        colors = np.array(colors)[keep_mask]

        # Default scatter style dictionary
        if scatter_style is None:
            stroke_width = self._stroke_width*0.375
            scatter_style = {"marker":"o",
                             "mstyle":{"stroke": "black",
                                       "stroke-width":stroke_width}}
            
        # Only plot if there there are nodes to plot. 
        if len(x) > 0:

            self._tree_ax.scatterplot(x,y,size=sizes,color=colors,**scatter_style)

            # If there was a property label, record that we plotted it for legend
            # construction.
            if property_label is not None and property_label.strip() != "":
                self._plotted_nodes[property_label] = (cm,cm_span,sm,sm_span,scatter_style)

        return self

    def draw_node_labels(self,
                         property_labels,
                         fmt_string=None,
                         position="right",
                         position_offset=None,
                         plot_ancestors=True,
                         plot_leaves=False,
                         plot_root=True,
                         text_style=None):
        """
        Draw labels on nodes on tree.

        Parameters
        ----------
        property_labels : str or list-like
            write values in property_labels to the node labels. if str, treat as
            a single property_label; if list-like, use each of the list entries
            as an individual property_label.
        fmt_string : str, optional
            format for labels. If not specified, spit out as {},{},... for each
            property_label entry. If specified, must have one of two formats.

            + exactly the same number of *unnamed* fields as (i.e. "{}") as the
              number of properties. For example, :code:`"{}|{}"` for
              :code:`property_labels = ["label_a","label_b"]`).
            + an arbitrary number of *named* fields (e.g. "{label_a}"). For
              example :code:`"{label_b}|{label_a}|{label_a}"` for
              :code:`property_labels = ["label_a","label_b"]`.

            Mixed format strings (i.e. :code:`"{}{label_a}"`) are not
            permitted, but all other python formatting should be supported.
        position : str, default="right"
            where to put label relative to node. allowed values are: top-left",
            "top", "top-right", "right", "bottom-right", "bottom", "bottom-left",
            and "left".
        position_offset : float, optional
            how far to displace the labels off the nodes, in px. If not specified,
            will be 0.75*node_size
        plot_ancestors : bool, default=True
            draw ancestral nodes
        plot_leaves : bool, default=False
            draw leaf nodes
        plot_root : bool, default=True
            draw the root node
        text_style : dict, optional
            dictionary specifying how to draw the text labels. toyplot Text
            documentation for description of allowable values.
        """

        # If passed in as a string, make into a list of strings
        if issubclass(type(property_labels),str):
            property_labels = [property_labels]

        # Check input bools
        plot_ancestors = check.check_bool(plot_ancestors,"plot_ancestors")
        plot_leaves = check.check_bool(plot_leaves,"plot_leaves")

        # Get position and position_offset
        position = str(position)
        if position_offset is not None:
            # Convert from px to tree coordinates.
            position_offset = check.check_float(position_offset,"position_offset")
            x_offset = np.abs(position_offset*self._x_px_to_tree)
            y_offset = np.abs(position_offset*self._y_px_to_tree)
        else:
            x_offset = np.abs(self.default_size*self._x_px_to_tree)
            y_offset = np.abs(self.default_size*self._y_px_to_tree)

        # small hack -- displace a bit more in y than x because it almost always
        # looks saner.
        y_offset = y_offset*1.25

        # Get displacement from node positions in x and y
        dx, dy = parse_position_string(position,x_offset,y_offset)

        # ----------------------------------------------------------------------
        # Deal with text style

        # Create text style dictionary
        if text_style is None:
            text_style = {}

        # Make sure text_style is sane
        if not issubclass(type(text_style),dict):
            err = "text_style should be a dictionary of css key/values\n"
            raise ValueError(err)

        # Load in defaults if user has not specified things like font-size
        try:
            text_style["font-size"]
        except KeyError:
            text_style["font-size"] = f"{self._font_size*0.75}px"

        try:
            text_style["text-anchor"]
        except KeyError:

            if np.isclose(dx,0):
                text_anchor = "middle"
            else:
                if dx < 0:
                    text_anchor = "end"
                else:
                    text_anchor = "start"

            text_style["text-anchor"] = text_anchor

        try:
            text_style["fill"]
        except KeyError:
            text_style["fill"] = color_to_css("#023E55")

        # ----------------------------------------------------------------------
        # Deal with fmt_string argument

        # Create default {label_a},{label_b},... style format string
        if fmt_string is None:
            fmt_string = ["{" + p + "}" for p in property_labels]
            fmt_string = ",".join(fmt_string)

        # Look for { } patterns in the fmt_string and try to match
        # them to property_label values. var_list holds variables to load in at
        # which position. cast_from holds the type from which this is being cast.
        var_list = []
        cast_from = []
        for m in re.finditer("{.*?}",fmt_string):

            # Get variable name from either {var_name} or {var_name:}
            v = m.string[m.start():m.end()]
            inside = v.strip("{").split("}")[0]

            # Get variable name
            fmt_bits = inside.split(":")
            var_name = fmt_bits[0]

            # See if this is numeric (has d or f in format)
            if len(fmt_bits) == 1:
                var_type = ""
            else:
                var_type = fmt_bits[1]

            number_fmt = re.search("[df]",var_type,flags=re.IGNORECASE)
            if number_fmt:
                match = number_fmt.group(0)
                if match == "d":
                    cast_from.append(int)
                else:
                    cast_from.append(float)
            else:
                cast_from.append(str)

            var_list.append(var_name)

            # Update fmt_string, wiping out var_name we just found
            fmt_string = re.sub(var_name,"",fmt_string)


        # no variable specified explicitly. if sane, will have exactly
        # len(property_labels) "" entries.
        if len(set(var_list)) == 1 and var_list[0] == "":
            if len(var_list) == len(property_labels):
                var_list = property_labels[:]

        # At this point, var_list should not have "" if fmt_string was sane.
        if "" in var_list:
            err = f"\nfmt_string '{fmt_string}' is ambiguous. It must either\n"
            err += "have explicitly indicated properties (i.e. {a_label:s}) or\n"
            err += "have *only* implicit labels (i.e. {}) with the same number\n"
            err += "of fields as requested property_labels.\n\n"
            raise ValueError(err)

        # Actually pull out x, y, and properties
        x, y, prop_dict = self._get_node_values(plot_ancestors,
                                                plot_leaves,
                                                plot_root,
                                                property_labels)

        # Offset in x and y
        x = x + dx
        y = y + dy

        good_mask = []

        # Apply formatting string to the property values extracted
        to_write = []
        to_zip = [prop_dict[var] for var in var_list]
        for v in zip(*to_zip):

            num_none = sum([value is None for value in v])
            if len(v) == num_none:
                good_mask.append(False)
                to_write.append(None)
                continue

            # Try to format
            try:
                new_string = fmt_string.format(*v)

            # If we hit a type error, problem with int or float. Turn float
            # into np.nan. Serious hack. Turn into '999999' and remove immediately
            # with re.sub. Ugly, but can't get custom formatter in here without
            # a huge refactor...
            except TypeError:

                # Go through each value in the vector
                v = list(v)
                for i in range(len(v)):

                    # see if cast_from works. If so, this is not problem value
                    try:
                        cast_from[i](v[i])

                    # If here, sub bad float with np.nan. sub bad int with
                    # 999999
                    except TypeError:
                        if issubclass(cast_from[i],float):
                            v[i] = np.nan
                        else:
                            v[i] = 999999

                # Format and get rid of wacky "999999" we added.
                new_string = fmt_string.format(*v)
                new_string = re.sub("999999","",new_string)

            good_mask.append(True)

            to_write.append(new_string)

        good_mask = np.array(good_mask,dtype=bool)
        to_write = np.array(to_write)

        # Actually write labels
        self._tree_ax.text(x[good_mask],
                           y[good_mask],
                           to_write[good_mask],
                           style=text_style)

        return self


    def draw_scale_bar(self,
                       bar_length=0.2,
                       units="subs/site"):
        """
        Draw a scale bar on the tree plot.

        Parameters
        ----------
        bar_length : float, default=0.3
            draw a scale bar that is bar_length fraction of the total tree
            length. Note: this may be ignored and the bar truncated if the
            bar would collide with the tree.
        units : str, default="subs/site"
            name to put below the scale bar.
        """

        bar_length = check.check_float(bar_length,
                                       "bar_length",
                                       minimum_allowed=0,
                                       maximum_allowed=1)

        # Get length of bar
        target_length = self._x_total_domain*bar_length
        if target_length > self._legend_width:
            target_length = self._legend_width

        # Round the branch length in a pretty fashion
        round_to = get_round_to(target_length)
        bar_length = np.round(target_length,round_to)
        if round_to == 0:
            bar_length = int(bar_length)



        bar_center = self._x_min + self._x_total_domain*0.15
        bar_y = 6

        # Plot scale bar
        bar_x_coord = (bar_center - bar_length/2,
                       bar_center + bar_length/2)
        bar_y_coord = (bar_y,bar_y)
        bar_tick_y = (bar_y-0.5,bar_y+0.5)

        self._legend_ax.plot(bar_x_coord,
                             bar_y_coord,
                             stroke_width=self._stroke_width*0.7,
                             color="black")
        self._legend_ax.plot((bar_x_coord[0],bar_x_coord[0]),
                             bar_tick_y,
                             stroke_width=self._stroke_width*0.7,
                             color="black")
        self._legend_ax.plot((bar_x_coord[1],bar_x_coord[1]),
                             bar_tick_y,
                             stroke_width=self._stroke_width*0.7,
                             color="black")

        # Label scale bar
        if round_to > 0:
            fmt = "{:." + str(round_to) + "f}"
        else:
            fmt = "{:d}"

        bar_label = f"{fmt.format(bar_length)} {units}"
        bar_label = bar_label.strip()
        bar_label_x = bar_center
        bar_label_y = bar_y - 2
        self._legend_ax.text(bar_label_x,
                             bar_label_y,
                             bar_label,
                             color="black",
                             style=self._draw_kwargs["tip_labels_style"])

        return self

    def draw_node_legend(self,
                         label_renamer=None,
                         max_label_len=25):
        """
        Add a node legend to the tree plot.

        Parameters
        ----------
        label_renamer : dict, optional
            dictionary where keys are property_labels for named series and
            values are what these should be called in the legend
        max_label_len : int, default=15
            truncate long property names to max_label_len
        """

        if label_renamer is None:
            label_renamer = {}
        if not issubclass(type(label_renamer),dict):
            err = "\nlabel_renamer '{label_renamer}' invalid. Should be a dict\n\n"
            raise ValueError(err)

        max_label_len = check.check_int(max_label_len,minimum_allowed=1)

        # Don't do anything if we don't have any named series
        if len(self._plotted_nodes) == 0:
            return self

        # Get label and node sizes for min and max of all nodes with plotted
        # properties.
        label_sizes = []
        spans = {}
        for node in self._plotted_nodes:

            # Get drawing information for the node series
            cm, cm_span, sm, sm_span, scatter_style = self._plotted_nodes[node]

            # Combine span for both cm and sm
            span = cm_span[:]
            span.extend(sm_span)
            span = list(set(span))
            span.sort()

            # If no span -- both don't vary -- set span to a single value
            if len(span) == 0:
                span = [1]

            # Prep label_renamer dictionary and get label length
            try:
                label = label_renamer[node]
            except KeyError:
                label_renamer[node] = node
                label = label_renamer[node]
            label_sizes.append(len(label))

            spans[node] = span

        # Overall legend node size will be the max seen across tree * 1.2 for
        # clarity.
        legend_node_size = self._font_size*1.2

        # Truncate label if necessary to accomodate max_label_len
        label_len = max(label_sizes)
        if label_len > max_label_len:
            label_len = max_label_len
        label_fmt = "{:>" + f"{label_len:d}" + "s}"

        # Get labels for all series, formatted correctly
        series_labels = []
        for node in self._plotted_nodes:
            label = label_renamer[node]
            if len(label) > label_len:
                series_labels.append(label_fmt.format(f"{label[:(label_len-3)]}..."))
            else:
                series_labels.append(label_fmt.format(label))

        # Now go through all plotted features
        for node_counter, node in enumerate(self._plotted_nodes):

            # Get node properties for plotting
            cm, _, sm, _, scatter_style = self._plotted_nodes[node]

            span = spans[node]

            # Construct generic marker_style that will be used for all markers
            marker_style = copy.deepcopy(scatter_style)
            marker_style["shape"] = marker_style["marker"]
            marker_style.pop("marker")
            marker_style["lstyle"] = {"font-size":f"{self._font_size*0.5}px"}
            marker_style["angle"] = 0

            # Figure out level at which to round value for marker label
            try:
                value = float(span[0])
                round_at = get_round_to(value)
            except (ValueError,TypeError):
                round_at = None

            # Create markers
            m_list = []
            for i in range(len(span)):

                # Format property value appropriately
                if round_at is None:
                    p = f"{span[i]}"
                elif round_at == 0:
                    p = f"{int(np.round(span[i],0)):d}"
                else:
                    prop_fmt = "{:." + f"{round_at}" + "f}"
                    p = prop_fmt.format(np.round(span[i],round_at))

                # Construct a style for this specific marker
                m_style = copy.deepcopy(marker_style)
                m_style["label"] = p
                m_style["mstyle"]["fill"] = cm(span[i])

                # Convert color to rgba vector
                rgba = toyplot.color.broadcast(m_style["mstyle"]["fill"],1)[0]
                saturation = (rgba[0] + rgba[1] + rgba[2])*rgba[3]
                if saturation < 1.5:
                    txt_color = "#ffffff"
                else:
                    txt_color = "#000000"

                m_style["lstyle"]["fill"] = txt_color
                m_style["size"] = legend_node_size

                m_list.append(toyplot.marker.Marker(**m_style).to_html())

            if round_at is None:
                sep = "  "
            else:
                sep = " &#8594; "

            marker_str = sep.join(m_list)
            marker_style = {"font-size":self._font_size,
                            "text-anchor":"start"}
            marker_x = self._x_min + self._x_total_domain*0.67
            self._legend_ax.text(marker_x,
                                 self._legend_y,
                                 marker_str,
                                 style=marker_style,
                                 color="black")

            label_style = {"font-size":self._font_size,
                            "text-anchor":"end"}
            label_str = f"{series_labels[node_counter]}:"
            label_x  = self._x_min + self._x_total_domain*0.65
            self._legend_ax.text(label_x,
                                 self._legend_y,
                                 label_str,
                                 style=label_style,
                                 color="black")

            self._legend_y -= 5

        return self

    def render(self,output_file):
        """
        Render a the tree out as a file.

        Parameters
        ----------
        output_file : str
            output file. should be .svg, .pdf, or .png
        """

        key = output_file[-3:].lower()
        render_dict = {"svg":toyplot.svg.render,
                       "pdf":toyplot.pdf.render,
                       "png":None}
        try:
            render_fcn = render_dict[key]
        except KeyError:
            print(f"Could not identify render type for file {output_file}.")
            print("Using pdf.",flush=True)
            render_fcn = render_dict["pdf"]

        # png rendering requires ghostscript, which requires a painful external
        # install chain (particularly on github workflows). Hack makes this a
        # fall back rather than something in the base import functionality.
        if key == "png":
            from toyplot import png
            render_fcn = toyplot.png.render

        render_fcn(self.canvas,output_file)

    def as_html(self):
        """
        Render the tree as an html string.

        Returns
        -------
        html : str
            tree as an html string
        """

        as_xml = toyplot.html.render(self.canvas)

        return xml.etree.ElementTree.tostring(as_xml, encoding='unicode')

    @property
    def canvas(self):
        """
        toyplot canvas holding tree.
        """
        return self._canvas

    @property
    def tree_ax(self):
        """
        toyplot axes holding tree.
        """

        return self._tree_ax

    @property
    def legend_ax(self):
        """
        toyplot axes holding legend.
        """

        return self._legend_ax

    @property
    def tree_mark(self):
        """
        tree object (a toyplot mark)
        """
        return self._tree_mark

    @property
    def tT(self):
        """
        toytree object used to generate the plot.
        """
        return self._tT

    @property
    def default_size(self):
        """
        Default size for nodes.
        """
        return self._stroke_width*6

    @property
    def plotted_properties(self):
        """
        Return list of named node series that have been plotted.
        """

        return list(self._plotted_nodes.keys())
