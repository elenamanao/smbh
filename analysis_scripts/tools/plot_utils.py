import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
# import seaborn as sns


# # Use this dictionary to have uniform colors between plots!
# source_color = {
#     'NGC1068': '#77B3D7',
#     'NGC4151': '#e28387',
#     'CGCG420-015': '#C01B3C',
#     # 'NGC 6240': '#785ef0',
#     # 'NGC 4388': 'tab:cyan',
#     'NGC7469': '#557687',
# }

# blues = ["#0C355D", "#1B5B99", "#789BC3"]
# azures = ["#557687", "#8CC0DA", "#BAD9E9"]
# greys = ["#42464F", "#8991A2", "#B7BCC6"]
# reds = ["#5F0819", "#C01B3C", "#e28387"]
# yellows = ["#87630D", "#DAA21C", "#EAC88A"]
# shades = ["#E4E986", "#8CDB5E", "#00C8CF"]
# all_colors = blues + azures + greys + reds + yellows
# palette_dark = all_colors[::3]
# palette = all_colors[1::3]
# palette_light = all_colors[2::3]

# pst_color = reds[1]
# ana_9y_color = yellows[1]
# ana_13y_color = blues[1]

# sns.set_palette(palette)

# marker_cycler = itertools.cycle(["o", "v", "^", "s", " D"])

nice_colors = colormaps["tab20b"]

colors = nice_colors(np.linspace(0, 1, 20))

blues = colors[0:4]
greens = colors[4:8]
yellows = colors[8:12]
reds = colors[12:16]
purples = colors[16:]


settings = {
    # Use LaTeX to write all text
    "text.usetex": True,
    # "pgf.preamble": r"\usepackage{amsmath}",
    # "backend": "pgf",
    # "pgf.texsystem": "pdflatex",
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
    # "svg.fonttype": "none",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.major.pad": 6,
    "ytick.major.pad": 6,
    "xtick.minor.bottom": False,
    "ytick.minor.left": False,
    "ytick.minor.right": False,
    # Custom axes
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    "axes.labelpad": 8,
    # Lines and markers
    "lines.linewidth": 3,
    'lines.markersize': 6,
    # Grids
    "grid.linestyle": ":",
    # Legend
    "legend.frameon": False,
    # Save figure
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
}

plt.rcParams.update(settings)


def set_size(width, fraction=0.8, subplots=(1, 1), ratio = None, height_ratio =[1]):
    """Set figure dimensions to avoid scaling in LaTeX.
    
    For PRX document I got the following widths:
    > 510.0pt.
    l.65 \showthe\textwidth

    > 246.0pt.
    l.66 \showthe\columnwidth

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predefined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    ratio: float or None
            Height / width ratio. If set to None, the golden ratio
            (5**0.5 - 1) / 2 is used.

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == "text":
        width_pt = 469.755
        settings = {
            "axes.labelsize": 10,
            "font.size": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "lines.linewidth": 2,
            "xtick.major.pad": 5,
            "ytick.major.pad": 5,
            "axes.labelpad": 8,
        }
        plt.rcParams.update(settings)
    elif width == "margin":
        width_pt = 144.54
        settings = {
            "axes.labelsize": 7,
            "font.size": 7,
            "legend.fontsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "lines.linewidth": 1.5,
            "xtick.major.pad": 3.5,
            "ytick.major.pad": 3.5,
            "axes.labelpad": 6,
        }
        plt.rcParams.update(settings)
    elif width == "prx_column":
        width_pt = 246.0
        settings = {
            "axes.labelsize": 10,
            "font.size": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "lines.linewidth": 2,
            "xtick.major.pad": 5,
            "ytick.major.pad": 5,
            "axes.labelpad": 8,
        }
        plt.rcParams.update(settings)
    elif width == "prx_wide":
        width_pt = 510.0
        settings = {
            "axes.labelsize": 10,
            "font.size": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "lines.linewidth": 2,
            "xtick.major.pad": 4,
            "ytick.major.pad": 4,
            "axes.labelpad": 8,
        }
        plt.rcParams.update(settings)


    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    if ratio == "equal":
        ratio = 1.0

    if ratio is None:
        # Golden ratio to set aesthetic figure height
        # https://disq.us/p/2940ij3
        ratio = (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * ratio * (subplots[0]*height_ratio[0] / (subplots[1]*np.sum(height_ratio)))

    return (fig_width_in, fig_height_in)
