"""
Class for drawing formatted phylogenetic trees using toytree.
"""

import topiary
from topiary.draw._core import construct_colormap, construct_sizemap
from topiary.draw._core import ete3_to_toytree
import topiary._private.check as check

import toytree
import toyplot

import numpy as np

import re, copy


def _get_round_to(value,total_requested=3):
    """
    Figure out a natural place to round a float for creating useful but pretty
    strings.

    Parameters
    ----------
    value : float
        value to round
    total_requested : int, default=3
        target number of digits for rounding. If value is 1234.23 and
        total_requested is 3, this would round to 1234. If total_requested
        is 5, this would give 1234.23.

    Returns
    -------
    round_to : int
        what number to round value to prior to string conversion
    """

    # Convert value to string, dropping sign
    value_str = f"{np.abs(value)}"

    # Deal with exponents (1e+50, for example)
    if re.search("e",value_str):

        # Split on "e"
        exp_split = value_str.split("e")
        base = exp_split[0]
        exponent = int(exp_split[1])

        # Process base
        base_split = base.split(".")
        base_whole = list(base_split[0])
        if len(base_split) == 2:
            base_decimal = list(base_split[1])
        else:
            base_decimal = []

        # Positive exponent
        if exponent > 0:
            decimal_out = ["0" for _ in range(int(exponent))]
            grab = min([len(decimal_out),len(base_decimal)])
            decimal_out[:grab] = base_decimal[:grab]
            value_str = "".join(base_whole) + "".join(decimal_out)

        # Negative exponent
        else:
            decimal_out = ["0" for _ in range(-exponent)]
            grab = min([len(decimal_out),len(base_decimal)])
            for i in range(grab):
                decimal_out[-(i+1)] = base_decimal[i]
            decimal_out[-(grab + len(base_whole)):-grab] = base_whole[:]

            value_str = "0." + "".join(decimal_out)

    # Split value on "."
    value_split = f"{value_str}".split(".")

    # If there is no decimal, round at 0
    if len(value_split) == 1:
        round_at = 0
    else:

        # Consider whole and decimal portions of the float
        whole = value_split[0]
        decimal = value_split[1]

        non_zero_num_whole = len(whole.lstrip("0"))

        round_at = 0
        num_taken = non_zero_num_whole
        for counter, d in enumerate(decimal):

            # If we've seen at least one non-zero digit and we've reached
            # requested total digits, break
            if num_taken >= total_requested:
                if non_zero_num_whole > 0 or round_at > 0:
                    break

            num_taken += 1
            if d != "0":
                round_at = (counter + 1)

    return round_at


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
    vertical_pixels_per_taxon : float, default=20
        number of pixels to assign to each taxon when calculating figure
        height
    aspect_ratio_parameters : tuple, default=(0.02,0.47)
        these parameters used to relate figure aspect ratio to number of
        taxa by the following equation:
        :code:`aspect_ratio_parameters[0]*num_taxa + aspect_ratio_parameters[1]`
    min_height : float, default=300
        minimum height for figure (pixels)
    **kwargs : dict
        pass any other keyword arguments directly to toytree.tree.draw
    """

    def __init__(self,
                 T,
                 name_dict=None,
                 font_size=15,
                 stroke_width=2,
                 vertical_pixels_per_taxon=20,
                 aspect_ratio_parameters=(0.02,0.47),
                 min_height=300,
                 **kwargs):
        """
        Initialize PrettyTree class.
        """

        # Read into an ete3 tree
        if issubclass(type(T),toytree.tree):
            self._tT = T.copy()
        else:
            # Convert to a toytree
            self._tT = ete3_to_toytree(topiary.io.read_tree(T))

        # Rename tree.name entries according to name_dict
        if name_dict is not None:
            for idx in self._tT.idx_dict:
                name = self._tT.idx_dict[idx].name
                try:
                    self._tT.idx_dict[idx].name = name_dict[name]
                except KeyError:
                    pass

        # Artwork parameters
        self._font_size = check.check_float(font_size,
                                            "font_size",
                                            minimum_allowed=0)
        self._stroke_width = check.check_float(stroke_width,
                                               "stroke_width",
                                               minimum_allowed=0)
        self._vertical_pixels_per_taxon = check.check_float(vertical_pixels_per_taxon,
                                                            "vertical_pixels_per_taxon",
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

        # If height not passed directly, figure out
        if "height" not in kwargs:
            num_nodes = self._tT.ntips
            self._height = num_nodes*self._vertical_pixels_per_taxon
            if self._height < min_height:
                self._height = min_height
            kwargs["height"] = f"{self._height}px"
        else:
            self._height = kwargs["height"]

        # If width not passed directly, figure out
        if "width" not in kwargs:
            ratio = aspect_ratio_parameters[0]*num_nodes + aspect_ratio_parameters[1]
            self._width = self._height/ratio
            kwargs["width"] = f"{self._width}px"
        else:
            self._width = kwargs["width"]

        # If padding not passsed directly, figure out
        if "padding" not in kwargs:
            kwargs["padding"] = 4*self._vertical_pixels_per_taxon

        # set stroke and stroke-width in edge_style
        if "edge_style" not in kwargs:
            kwargs["edge_style"] = {}
        kwargs["edge_style"]["stroke-width"] = self._stroke_width
        if "stroke" not in kwargs["edge_style"]:
            kwargs["edge_style"]["stroke"] = "black"

        # set font-size in tip_labels_style
        if "tip_labels_style" not in kwargs:
            kwargs["tip_labels_style"] = {}
        kwargs["tip_labels_style"]["font-size"] = f"{self._font_size}px"

        # Get final kwargs
        self._draw_kwargs = kwargs

        # Draw the first draft of the tree
        self._canvas, self._axes, self._mark = self._tT.draw(**self._draw_kwargs)

        # Dictionary to hold named nodes added
        self._named_node_series = {}

        # Get drawing attributes we'll need
        self._get_pixel_aspect_ratio()
        self._get_dimensions()
        self._get_legend_location()

    def _get_dimensions(self):
        """
        Get the plot/tree dimensions.
        """

        # Size in tree coordinates
        self._min_x, self._max_x = self._mark.domain("x")
        self._min_y, self._max_y = self._mark.domain("y")

        node_coord = self._tT.get_node_coordinates()

        # Get indexes for furthest left and furthest right nodes
        min_x_idx = np.argmin(node_coord[:,0])
        max_x_idx = np.argmax(node_coord[:,0])

        # Get all edges
        edges = self._tT.get_edges()

        # Get total length of tree
        self._total_length = self._max_x - self._min_x

        # Separation between tips
        self._y_step = -self._font_size*1.333*self._y_px_to_tree #(self._max_y - self._min_y)/self._tT.ntips


    def _get_pixel_aspect_ratio(self):
        """
        Get the aspect ratio of each pixel (toyplot domain) in terms of the
        tree coordinate (toyplot range).
        """

        # Draw nodes of size zero and render to allow domain (pixel) calcs
        # to work
        node_coord = self._tT.get_node_coordinates()
        self._axes.scatterplot(node_coord[:,0],node_coord[:,1],size=0)
        _ = toyplot.html.render(self.canvas)

        # Get map between domain and range in x
        x_domain = np.array(self._mark.domain("x"))
        x_range = self._axes.project("x",x_domain)
        x_domain_to_range = (x_domain[1] - x_domain[0])/(x_range[1] - x_range[0])
        self._x_px_to_tree = x_domain_to_range
        self._x_tree_to_px = 1/x_domain_to_range

        # Get map between domain and range in y
        y_domain = np.array(self._mark.domain("y"))
        y_range = self._axes.project("y",y_domain)
        y_domain_to_range = (y_domain[1] - y_domain[0])/(y_range[1] - y_range[0])
        self._y_px_to_tree = y_domain_to_range
        self._y_tree_to_px = 1/y_domain_to_range

        self._pixel_aspect = np.abs(x_domain_to_range/y_domain_to_range)


    def _get_legend_location(self):
        """
        Get the location to place legend based on tree topology.
        """

        # get node coordinates
        node_coord = self._tT.get_node_coordinates()

        # Create a dictionary chaining edges backwards
        edges = self._tT.get_edges()
        edge_dict = dict(zip(edges[:,1],edges[:,0]))

        # Node at bottom of graph (lowest y) -- will be a leaf
        start_idx = np.argmin(node_coord[:,1])

        # Node at left of graph (lowest x) -- also ancestor
        end_idx = np.argmin(node_coord[:,0])

        # Bottom left corner of the graph
        corner = np.array([self._min_x,self._min_y])

        best_aspect_node = None
        best_aspect_diff = np.inf

        # Walk through the chain of nodes starting at min_y, finding node that
        # draws the closest thing to a square with the bottom left corner
        node_idx = start_idx
        while node_idx != end_idx:

            # If y displacement is 0, make very low...
            diff = np.abs(node_coord[node_idx,:] - corner)
            if diff[1] == 0:
                diff[1] = 1e-6

            aspect = (diff[0]/diff[1])/self._pixel_aspect
            aspect_diff = np.abs(1 - aspect)
            if aspect_diff < best_aspect_diff:
                best_aspect_diff = aspect_diff
                best_aspect_node = node_coord[node_idx,:]

            # Move back a node. If KeyError, tree is in strange configuration
            # and cannot be parsed visually this way
            try:
                node_idx = edge_dict[node_idx]
            except KeyError:
                break

        # If we did not find best_aspect_node
        if best_aspect_node is None:
            x = (self._max_x - self._min_x)/8 + self._min_x
            y = (self._max_y - self._min_y)/8 + self._min_y
            best_aspect_node = [x,y]

        box_width = np.abs(corner[0] - best_aspect_node[0])
        box_height = np.abs(corner[1] - best_aspect_node[1])

        self._legend_x = (corner[0] + best_aspect_node[0])/2
        self._legend_y = corner[1] + self._y_step
        self._legend_width = box_width

    def _get_node_values(self,get_ancestors,get_leaves,get_root=True,property_label=None):
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
        property_label : str or list-like
            list of property labels to grab, in order. can take a single string
            as equivalent to a length-1 list.

        Returns
        -------
        x : numpy.ndarray
            x coordinates of nodes
        y : numpy.ndarray
            y coordinates of nodes
        all_props : dict or None
            dictionary keying property_label to an np.ndarray of that property
            value. if property_label is None, return None.
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
        if property_label is not None:

            # If passed in as a single string, put in a list to allow iteration
            # over single label
            if issubclass(type(property_label),str):
                property_label = [property_label]

            # Iterate over all labels
            for p_label in property_label:

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

                # Figure out how to round value
                if len(prop) > 0:
                    prop = np.array(prop)
                    if np.issubdtype(type(prop[0]), np.floating):
                        min_round_to = _get_round_to(np.min(prop))
                        max_round_to = _get_round_to(np.max(prop))
                        round_to = np.max([min_round_to,max_round_to])

                        if round_to == 0:
                            prop = np.array(np.round(prop,0),dtype=int)
                        else:
                            prop = np.round(prop,round_to)

                all_props[p_label] = prop

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
        prop = np.array(prop)
        x = x[good_mask]
        y = y[good_mask]
        prop = prop[good_mask]

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

        cm, cm_span = construct_colormap(color,prop,prop_span,palette)

        # Construct size map
        if size is None:
            size = self.default_size

        sm, sm_span = construct_sizemap(size,prop,prop_span)

        # Construct size and color lists for each node
        sizes = []
        colors = []
        for p in prop:
            sizes.append(sm(p))
            colors.append(cm(p))

        # Default scatter style dictionary
        if scatter_style is None:
            stroke_width = self._stroke_width*0.375
            scatter_style = {"marker":"o",
                             "mstyle":{"stroke": "black",
                                       "stroke-width":stroke_width}}

        self._axes.scatterplot(x,y,size=sizes,color=colors,**scatter_style)

        # If there was a property label, record that we plotted it for legend
        # construction.
        if property_label is not None and property_label.strip() != "":
            self._named_node_series[property_label] = (cm,cm_span,sm,sm_span,scatter_style)

        return self

    def draw_node_labels(self,
                         property_labels,
                         fmt_string=None,
                         plot_ancestors=True,
                         plot_leaves=False,
                         plot_root=True,
                         text_style=None,
                         node_size=None,
                         draw_below=False):
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
            + an arbitrary number of *named* fields (e.g. "{label_a"). For
              example :code:`"{label_b}|{label_a}|{label_a}"` for
              :code:`property_labels = ["label_a","label_b"]`.
            Mixed format strings (i.e. :code:`"{}{label_a}"`) are not
            permitted, but all other python formatting should be supported.
        plot_ancestors : bool, default=True
            draw ancestral nodes
        plot_leaves : bool, default=False
            draw leaf nodes
        plot_root : bool, default=True
            draw the root node
        text_style : dict, optional
            dictionary specifying how to draw the text labels. toyplot Text
            documentation for description of allowable values.
        node_size : float, optional
            draw label assuming a node of size node_size. if None, use the
            default node size (stroke_width*6)
        draw_below : bool, default=False
            draw label below the branch. otherwise, draw above.
        """

        # Make property_label into a list-like if single string
        if issubclass(type(property_labels),str):
            property_labels = [property_labels]

        plot_ancestors = check.check_bool(plot_ancestors,"plot_ancestors")
        plot_leaves = check.check_bool(plot_leaves,"plot_leaves")
        draw_below = check.check_bool(draw_below,"draw_below")

        x, y, prop_dict = self._get_node_values(plot_ancestors,
                                                plot_leaves,
                                                plot_root,
                                                property_labels)

        if node_size is None:
            node_size = self._stroke_width*6
        node_size = check.check_float(node_size,"node_size",minimum_allowed=0)

        # Offset x
        x = x - self._x_px_to_tree*node_size*0.5

        # Offset y
        if draw_below:
            y = y - self._y_step*0.35
        else:
            y = y + self._y_step*0.35

        # Create text style dictionary
        if text_style is None:
            text_style = {"font-size":f"{self._font_size*0.75}px",
                          "text-anchor":"end", # right-justify
                          "fill":"black"}

        # Create default {label_a},{label_b},... style format string
        if fmt_string is None:
            fmt_string = ["{" + p + "}" for p in property_labels]
            fmt_string = ",".join(fmt_string)

        # Look for {xxx} patterns in the fmt_string and try to match
        # them to property_label values.
        var_list = []
        for m in re.finditer("{.*?}",fmt_string):

            # Get variable name from either {var_name} or {var_name:}
            v = m.string[m.start():m.end()]
            var_name = v.strip("{").split("}")[0].split(":")[0]
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

        # Apply formatting string to the property values extracted
        to_write = []
        to_zip = [prop_dict[var] for var in var_list]
        for v in zip(*to_zip):
            to_write.append(fmt_string.format(*v))

        # Actually write labels
        self._axes.text(x,y,to_write,style=text_style)

        return self


    def draw_scale_bar(self,
                       bar_length=0.3,
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
        target_length = self._total_length*bar_length
        if target_length > self._legend_width:
            target_length = self._legend_width

        # Round the branch length in a pretty fashion
        round_to = _get_round_to(target_length)
        bar_length = np.round(target_length,round_to)
        if round_to == 0:
            bar_length = int(bar_length)

        # Plot scale bar
        bar_x_coord = (self._legend_x - bar_length/2,self._legend_x + bar_length/2)
        bar_y_coord = (self._legend_y,self._legend_y)
        bar_tick_y = (self._legend_y-self._y_step*0.2,self._legend_y+self._y_step*0.2)
        self._axes.plot(bar_x_coord,
                        bar_y_coord,
                        stroke_width=self._stroke_width*0.5,
                        color="black")
        self._axes.plot((bar_x_coord[0],bar_x_coord[0]),
                        bar_tick_y,
                        stroke_width=self._stroke_width*0.5,
                        color="black")
        self._axes.plot((bar_x_coord[1],bar_x_coord[1]),
                        bar_tick_y,
                        stroke_width=self._stroke_width*0.5,
                        color="black")

        # Label scale bar
        if round_to > 0:
            fmt = "{:." + str(round_to) + "f}"
        else:
            fmt = "{:d}"

        bar_label = f"{fmt.format(bar_length)} {units}"
        bar_label = bar_label.strip()
        bar_label_x = self._legend_x
        bar_label_y = self._legend_y - 0.7*self._y_step
        self._axes.text(bar_label_x,
                        bar_label_y,
                        bar_label,
                        color="black",
                        style=self._draw_kwargs["tip_labels_style"])

        # Move legend_y down for next legend element.
        self._legend_y -= 2*self._y_step

        return self

    def draw_node_legend(self,
                         label_renamer=None,
                         max_label_len=15):
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

        # Get label and node sizes for min and max of all nodes with plotted
        # properties.
        label_sizes = []
        spans = {}
        for node in self._named_node_series:

            # Get drawing information for the node series
            cm, cm_span, sm, sm_span, scatter_style = self._named_node_series[node]

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

        # Overall legend node size will be the max seen across tree * 1.5 for
        # clarity.
        legend_node_size = self._font_size*1.2

        # Truncate label if necessary to accomodate max_label_len
        label_len = max(label_sizes)
        if label_len > max_label_len:
            label_len = max_label_len
        label_fmt = "{:>" + f"{label_len:d}" + "s}"

        # Get labels for all series, formatted correctly
        series_labels = []
        for node in self._named_node_series:
            label = label_renamer[node]
            if len(label) > label_len:
                series_labels.append(label_fmt.format(f"{label[:(label_len-3)]}..."))
            else:
                series_labels.append(label_fmt.format(label))

        # Now go through all plotted features
        for node_counter, node in enumerate(self._named_node_series):

            # Get node properties for plotting
            cm, _, sm, _, scatter_style = self._named_node_series[node]

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
                round_at = _get_round_to(value)
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

            m_str = sep.join(m_list)
            txt_str = f"{series_labels[node_counter]}: {m_str}"

            style = {"font-size":self._font_size}
            self._axes.text(self._legend_x,
                            self._legend_y,
                            txt_str,
                            style=style,
                            color="black")

            self._legend_y -= 1.5*self._y_step

        return self


    @property
    def canvas(self):
        """
        toyplot canvas holding tree.
        """
        return self._canvas

    @property
    def axes(self):
        """
        toyplot axes holding tree.
        """

        return self._axes

    @property
    def mark(self):
        """
        tree object (a toyplot mark)
        """
        return self._mark

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
