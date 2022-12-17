"""
Functions for generating html elements used in the report.
"""
import topiary
from topiary.draw.core import construct_colormap

import pandas as pd
import toyplot
import numpy as np

import xml
import os
import shutil
import copy

def create_output_directory(output_directory,overwrite=False):
    """
    Create an output directory for a topiary report, copying in assets.
    
    Parameters
    ----------
    output_directory : str
        directory to create
    overwrite : bool, default=False
        whether or not to overwrite an existing directory
    """

    # Create output directory
    if os.path.isdir(output_directory):
        if overwrite:
            shutil.rmtree(output_directory)
        else:
            err = f"output_directory '{output_directory}' already exists. Either\n"
            err += "choose a new output directory or set overwrite = True.\n\n"
            raise FileExistsError(err)
    os.mkdir(output_directory)

    # Copy in asset directory
    asset_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "assets")
    shutil.copytree(asset_dir,os.path.join(output_directory,".assets"))


def create_main_html(description="",
                     title="",
                     custom_css=".assets/topiary.css",
                     bs_css="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css",
                     bs_css_key="sha384-1BmE4kWBq78iYhFldvKuhfTAU6auU8tT94WrHftjDbrCEXSU1oBoqyl2QvZ6jIW3",
                     bs_js="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js",
                     bs_js_key="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p"):
    """
    Create the core elements of an html file. Top starts with <!doctype html>
    and ends with <body>. Bottom starts with <script></body> and ends with </html>.
    top + bottom will create a valid (empty) html file.

    Parameters
    ----------
    description : str
        text for meta description element
    title : str
        text for meta title element
    custom_css : str, optional
        path to a custom css file
    bs_css : str
        cdn url to bootstrap css file 
    bs_css_key : str
        key used to access cdn
    bs_js : str
        cdn url to bootstrap min js file
    bs_js_key : str
        key used to access cdn

    Returns
    -------
    top : str
        top half of a valid html file (ends with <body>)
    bottom : str
        bottom half of a valid html file (starts with <script></body>)
    """


    out = ["<!doctype html>"]
    out.append("<html lang=\"en\">")
    out.append("<head>")
    out.append( "  <meta charset=\"utf-8\">")
    out.append( "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">")
    out.append(f"  <meta name=\"description\" content=\"{description}\">")
    out.append(f"  <meta name=\"generator\" content=\"topiary {topiary.__version__}\">")
    out.append(f"  <title>{title}</title>")
    out.append(f"  <link href=\"{bs_css}\" rel=\"stylesheet\" integrity=\"{bs_css_key}\" crossorigin=\"anonymous\">")
    out.append(f"  <link rel=\"shortcut icon\" type=\"image/x-icon\" href=\".assets/favicon.ico\">")

    if custom_css is not None:
        out.append(f"  <link href=\"{custom_css}\" rel=\"stylesheet\" type=\"text/css\">")
    
    out.append("</head>")
    out.append("")
    out.append("<body>")
    out.append("")

    top = "\n".join(out)

    out = [""]
    out.append(f"<script src=\"{bs_js}\" integrity=\"{bs_js_key}\" crossorigin=\"anonymous\"></script>")
    out.append("</body>")
    out.append("</html>")
    bottom = "\n".join(out)

    return top, bottom


def df_to_table(df,add_header=True,show_row_numbers=True,float_fmt="{}",int_fmt="{}"):
    """
    Generate an html table given a pandas dataframe.
    
    Parameters
    ----------
    df : pandas.DataFrame
        pandas dataframe to convert to an html table
    add_header : bool, default=True
        whether or not to have the header row on the table.
    show_row_numbers : bool, default=True
        whether or not to show row numbers on the table. 
    float_fmt : str, default="{}"
        format all float numbers using this format string
    int_fmt : str, default="{}"
            format all int numbers using this format string

    Returns
    -------
    out : str
        html holding table
    """

    if not issubclass(type(df),pd.DataFrame):
        err = "df should be a pandas DataFrame\n"
        raise ValueError(err)

    columns = list(df.columns)

    header = ["<tr>"]

    if show_row_numbers:
        header.append('<th scope="col">#</th>')

    for col in columns:
        header.append('<th scope="col">')
        header.append(col)
        header.append("</th>")
    header.append("</tr>")

    header = "".join(header)

    contents = []
    for i in range(len(df.index)):
        contents.append("<tr>")

        if show_row_numbers:
            contents.append(f"<th scope='row'>{i+1}</th>")

        for col in columns:
            
            value = df.loc[df.index[i],col]
            if np.issubdtype(type(value),int):
                value = int_fmt.format(value)
            elif np.issubdtype(type(value),float):
                value = float_fmt.format(value)
            else:
                value = "{}".format(value)

            contents.append(f"<td>{value}</td>")
        contents.append("</tr>")
    
    contents = "".join(contents)

    final_out = []
    final_out.append('<div class="table-responsive">')
    final_out.append('<table class="table table-striped table-sm table-hover">')

    if add_header:
        final_out.append(f"<thead>{header}</thead>")

    final_out.append(f"<tbody>{contents}</tbody>")
    final_out.append("</table>")
    final_out.append("</div>")

    final_out = "".join(final_out)
    return f"<div class=\"text-center align-middle\">{final_out}</div>"

def canvas_to_html(toyplot_canvas):
    """
    Render a toyplot canvas as html.
    
    Parameters
    ----------
    toyplot_canvas : toyplot.canvas
        object to render as html
    
    Returns
    -------
    out : str
        html holding tree
    """

    as_xml = toyplot.html.render(toyplot_canvas).findall('svg')
    
    out = []
    for s in as_xml:
        out.append(xml.etree.ElementTree.tostring(s, encoding='unicode'))
    
    html = "".join(out)

    return f"<div class=\"text-center\">{html}</div>"

def sequence_box(text,
                 color="#000000",
                 prop_value=None,
                 prop_span=None,
                 palette=None):
    """
    Construct a sequence box (inside a bootstrap card div) where each letter in 
    the text is potentially colored by prop_value.

    Parameters
    ----------

    text : list-like
        sequence to write out
    color : str or tuple or dict
        set text color. If a single value, color all letters that color. If
        list-like and length 2, treat as colors for minimum and maximum of a
        color gradient.  If dict, map property keys to color values. Colors
        can be RGBA tuples, named colors, or hexadecimal strings. See the
        toyplot documentation for allowed values.
    prop_value : list-like
        values of properties over which to make the color map. 
    prop_span : tuple, optional
        set min/max values for property color/size calculation. First element
        is min, second is max. Only applicable if color is quantitative gradient.
    palette : toyplot.color.Palette, optional
        custom toyplot palette. if specified, use this palette rather than
        color to define color scheme for a color gradient. requires the
        property value be a float.

    Returns
    -------
    out : str
        html with bootstrap card div holding sequence
    """

    out = []

    start, _ = create_element("div",{"class":["card","seq-box"]})
    out.append(start)

    start, _ = create_element("div",{"class":["card-body",
                                              "seq-box",
                                              "font-monospace"]})
    out.append(start)
    
    # If both of these conditions are met, we need to contruct a text string 
    # with a set of spans
    if prop_value is not None and not issubclass(type(color),str):

        if len(text) != len(prop_value):
            err= "text and prop_value must be the same length.\n"
            raise ValueError(err)

        cm, _ = construct_colormap(color=color,
                                   prop=prop_value,
                                   prop_span=prop_span,
                                   palette=palette)

        new_text = []
        for i, t in enumerate(list(text)):  
            this_color = cm(prop_value[i])
            new_text.append(f"<span style=\"color:{this_color}\">{t}</span>")
        
        text = "".join(new_text)

    else:
        text = "".join(text)
        text = f"<span style=\"color:{color}\">{text}</span>"

    out.append(text)
        
    out.append("</div></div>")

    return "".join(out)


def create_card(card_title=None,card_contents=None,title_tag="h6",match_height=True):
    """
    Create a bootstrap card with title and contents. 

    Parameters
    ----------
    card_title : str, optional
        text for card title. if None, do not create a title element
    card_contents : str, optional
        text for card contents. if None, do not add card contents
    title_tag : str, default="h6"
        text for title (h1, p, etc.)
    match_height : bool, default = True
        whether or not to match the height of this card to adjacent cards

    Returns
    -------
    out : str
        html card element
    """

    out = []

    if match_height:
        out.append("<div class=\"card h-100\">")
    else:
        out.append("<div class=\"card\">")

    out.append("<div class=\"card-body\">")
    if card_title is not None:
        out.append(f"<{title_tag} class=\"card-title\">{card_title}</{title_tag}>")
    if card_contents is not None:
        out.append(f"<div class=\"align-middle\">{card_contents}</div>")
    out.append("</div></div>")

    return "".join(out)

def create_element(element,attributes=None):
    """
    Create a generic html element. 
    
    Parameters
    ----------
    element : str
        element type (e.g., h1, button, etc.)
    attributes : dict
        element attributes as key/value pairs (e.g. {"id":"element_id", etc.)

    Returns
    -------
    start : str
        opening tag (e.g., <element attr1="1" attr2="2">)
    end : str
        closing tag (e.g., </element>)
    """

    out = [f"<{element}"]

    if attributes is not None:
        attributes = copy.deepcopy(attributes)

        for a in attributes:
            if hasattr(attributes[a],"__iter__"):
                if issubclass(type(attributes[a]),str):
                    attributes[a] = [attributes[a]]
            
            if len(attributes[a]) < 1:
                continue

            out.append(f" {a}=\"")
            out.append(" ".join([f"{v}" for v in attributes[a]]))
            out.append("\"")
    
    out.append(">")

    return "".join(out), f"</{element}>"

def create_icon_row(files_to_link,descriptions):
    """
    Create a row of icons.
    
    Parameters
    ----------
    files_to_link : list
        list of files to create link icons for
    descriptions : list
        list of descriptions to apply as tooltips to icons

    Returns
    -------
    out : str
        html element encoding icon row
    """

    ext_files = {"csv":".assets/csv_icon.svg",
                 "pdf":".assets/pdf_icon.svg",
                 "txt":".assets/txt_icon.svg"}

    out = []
    for i, f in enumerate(files_to_link):

        try:
            icon = ext_files[f[-3:]]
        except KeyError:
            icon = ext_files["txt"]

        s, e = create_element("a",{"class":["btn","btn-outline-dark"],
                                    "data-toggle":"tooltip",
                                    "data-html":"true",
                                    "title":descriptions[i],
                                    "href":f})

        out.append(s)
        out.append(f"<img src=\"{icon}\" class=\"img-fluid\" width=\"35px\" height=\"35px\" />")
        out.append(e)
    
    return "".join(out)

def create_row(elements):
    """
    Create a row of columns with this format

    <div class="text-center">
        <div class="row">
            <div class="col">
                e0
            </div>
            <div class="col">
                e1
            </div>
            ...
        </div>
    </div>

    Parameters
    ----------
    elements : list
        list of strings to put in as columns

    Returns
    -------
    html : str
        html result
    """

    con_s, con_e = create_element("div",{"class":["text-center","align-middle"]})
    rs, re = create_element("div",{"class":["row"]})
    cs, ce = create_element("div",{"class":["col"]})

    this_out = [f"{con_s}{rs}"]
    for e in elements:
        this_out.append(f"{cs}{e}{ce}")
    this_out.append(f"{re}{con_e}")
    
    return "".join(this_out)
    
        
        
        