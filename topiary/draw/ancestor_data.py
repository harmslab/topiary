"""
Create a summary plot for an ancestor.
"""

import numpy as np

from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec

# -----------------------------------------------------------------------------
# Configure plotting
# -----------------------------------------------------------------------------

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def _draw_histogram(values,ax,bin_size=0.05,color="gray"):
    """
    Draw a histogram sideways next to main plot.

    Parameters
    ----------
    values : numpy.array (float)
        array to bin
    ax : matplotlib.Axes
        axis on which to make plot
    bin_size : float
        width of bins (from 0 to 1.0)
    color : str
        color to make histogram bars

    Returns
    -------
    max_counts : int
        maximum number of counts (for constructing xlim later)
    """

    # Create histogram
    counts, bins = np.histogram(values,
                                bins=np.arange(0,1 + bin_size,bin_size))

    # Draw bars for histogram
    for i in range(len(counts)):
        rect = patches.Rectangle((0,bins[i]),
                                 width=counts[i],
                                 height=(bins[i+1]-bins[i]),
                                 linewidth=1,
                                 edgecolor='black',
                                 facecolor=color,
                                 alpha=0.5)
        ax.add_patch(rect)

    return np.max(counts)

def plot_ancestor_data(df_anc,
                       alt_anc_pp=0.25,
                       width_ratio=5,
                       anc_name=None,
                       anc_data_string=None):
    """
    Create a summary plot for an ancestor.

    Parameters
    ----------
    df_anc : pandas.DataFrame
        ancestral data frame
    alt_anc_pp : float, default=0.25
        cutoff (inclusive) for identifying plausible alternate states
    width_ratio : float, default=5
        width ratio for plot
    anc_name : str, optional
        name of ancestor (title on graph and output file {anc_name}.pdf)
    anc_data_string : str, optional
        data to dump in subtitle

    Returns
    -------
    fig : matplotlib.Figure
        matplotlib figure instance
    axes : numpy.array
        array of matplotlib axes (0 is main plot; 1 is histogram)
    """

    # Data frames containing unambiguous gaps and other sites
    df_gap = df_anc.loc[df_anc.site_type == "gap",:]
    df_nogap = df_anc.loc[df_anc.site_type != "gap",:]

    # create a figure
    fig = plt.figure()
    fig.set_figheight(4)
    fig.set_figwidth(10)

    # create grid for main and histogram subplots
    spec = gridspec.GridSpec(ncols=2, nrows=1,
                             width_ratios=[width_ratio, 1],
                             wspace=0.01)
    # Greate actual subplots
    ax = [fig.add_subplot(spec[0])]
    ax.append(fig.add_subplot(spec[1],sharey=ax[0]))

    # Plot gaps on main figure. Draw boxes for contiguous gap regions, Code
    # below is too clever by half, but returns contiguous blocks of gaps
    sites = np.array(df_gap.site,dtype=np.uint)
    contiguous = np.split(sites,np.where(np.diff(sites) != 1)[0]+1)

    # Iterate over contiguous blocks
    for c in contiguous:

        # Skip contigous gaps with no length.
        if len(c) == 0:
            continue

        # Find gap width.  Minimum gap width is 1.
        width = c[-1] - c[0]
        if width == 0:
            width = 1

        # Draw rectangle for gap block
        rect = patches.Rectangle((c[0],0),width,1,
                                 linewidth=1,
                                 edgecolor='lightgray',
                                 facecolor='lightgray')
        ax[0].add_patch(rect)

    # Create list of ambiguous gaps and draw veritical purple lines at these
    # positions.
    ambig_df = df_anc.loc[df_anc.site_type == "possible gap",:]
    for i in range(len(ambig_df)):
        row = ambig_df.iloc[i]
        ax[0].plot((row.site,row.site),(0,1.05),"--",lw=1,color="purple")

    # Plot ML and alt pp points
    ax[0].plot(df_nogap.site,df_nogap.ml_pp,".",color="black",markersize=8)
    ax[0].plot(df_nogap.site,df_nogap.alt_pp,".",color="red",markersize=8)

    # Plot ML and alt pp lines. Only draw lines over contiguous stretches.
    sites = np.array(df_nogap.site,dtype=np.uint)
    contiguous = np.split(sites,np.where(np.diff(sites) != 1)[0]+1)
    for c in contiguous:
        ax[0].plot(df_anc.site.iloc[c],df_anc.ml_pp.iloc[c],color="black",lw=2)
        ax[0].plot(df_anc.site.iloc[c],df_anc.alt_pp.iloc[c],color="red",lw=2)

    # Plot alt-all cutoff line
    ax[0].plot((np.min(df_anc.site),np.max(df_anc.site)),
               (alt_anc_pp,alt_anc_pp),"--",color="gray")

    # Draw histograms for ml and alt pp in right plot
    max_ml_counts = _draw_histogram(df_nogap.ml_pp,ax[1],color="gray")
    max_alt_counts = _draw_histogram(df_nogap.alt_pp,ax[1],color="red")

    # Figure out biggest value seen in histogram plot
    hist_x_max = np.max((max_ml_counts,max_alt_counts))

    # Plot alt-all cutoff line on histogram
    ax[1].plot((0,hist_x_max),
               (alt_anc_pp,alt_anc_pp),"--",color="gray")

    # Clean up axes for main plot
    ax[0].set_xlabel("alignment site")
    ax[0].set_ylabel("posterior probability")
    ax[0].set_xlim(np.min(df_anc.site),np.max(df_anc.site))
    ax[0].spines["right"].set_visible(False)
    ax[0].spines["top"].set_visible(False)
    ax[0].set_ylim(-0.05,1.1)

    # Clean up axes for histogram plot
    ax[1].set_xlim(0,hist_x_max*1.1)
    ax[1].spines["right"].set_visible(False)
    ax[1].spines["top"].set_visible(False)
    ax[1].spines["bottom"].set_visible(True)
    ax[1].spines["left"].set_visible(False)
    ax[1].axis('off')

    # Main plot title
    if anc_name is not None:
        fig.suptitle(f"{anc_name}")

    # This bit of wackiness finds where on main plot the midpoint for the
    # whole graph is.
    main_plot_x_length = np.max(df_anc.site) - np.min(df_anc.site)
    total_plot_length = ((1 + width_ratio)/(width_ratio))*main_plot_x_length
    mid_point = total_plot_length/2.0

    # Plot the subtitle on the main plot
    if anc_data_string is not None:
        ax[0].text(mid_point,1.08,anc_data_string,ha='center')

    # Save out figure
    if anc_name is not None:
        fig.savefig(f"{anc_name}.pdf",bbox_inches="tight")

    return fig, ax
