#from https://jwalton.info/Embed-Publication-Matplotlib-Latex/

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

import numpy as np
import scipy
import scipy.interpolate

import string

def get_figure_size(width, fraction=1, ratio=(5**.5 - 1) / 2, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 426.79135
    elif width == 'beamer':
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

def to_roman(num):
	res = ""
	table = [(1000, "m"), (900, "cm"), (500, "d"), (400, "cd"), (100, "c"), (90, "xc"), (50, "l"), (40, "xl"), (10, "x"), (9, "ix"), (5, "v"), (4, "iv"), (1, "i")]
	for cap, roman in table:
		d, m = divmod(num, cap)
		res += roman * d
		num = m
	return res

def add_label_ext_2d(ax, rowId, colId):
	label = list(string.ascii_lowercase)
	ax.set_title(label[rowId] + "." + to_roman(colId) + ")", loc='left', fontsize=8)

def add_label_ext(ax, labelId):
	label = list(string.ascii_lowercase)
	ax.set_title(label[labelId] + ")", loc='left', fontsize=8)

def add_label_int_2d(fig, ax, rowId, colId, facecolor=0.7):
	label = list(string.ascii_lowercase)
	
	trans = mtransforms.ScaledTranslation(2/72, -2/72, fig.dpi_scale_trans)
	ax.text(0.0, 1.0, label[rowId] + "." + to_roman(colId) + ")", transform=ax.transAxes + trans,
            fontsize=8, verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor=str(facecolor), edgecolor='none', pad=1.0))

def add_label_int(fig, ax, labelId, facecolor=0.7):
	label = list(string.ascii_lowercase)
	
	trans = mtransforms.ScaledTranslation(2/72, -2/72, fig.dpi_scale_trans)
	ax.text(0.0, 1.0, label[labelId] + ")", transform=ax.transAxes + trans,
            fontsize=8, verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor=str(facecolor), edgecolor='none', pad=1.0))

def interp2d(x, y, f, N, interp_type='cubic'):
	nx, ny = N
	
	xnew = np.linspace(x[0], x[-1], nx)
	ynew = np.linspace(y[0], y[-1], ny)
	
	interp_f = scipy.interpolate.interp2d(x, y, f, kind=interp_type)
	return xnew, ynew, interp_f(xnew, ynew)

# plt.style.use('seaborn')

tex_fonts = {
	# Use LaTeX to write all text
	"text.usetex": True,
	"font.family": "serif",
	# Use 10pt font in plots, to match 10pt font in document
	"axes.labelsize": 10,
	"font.size": 10,
	# Make the legend/label fonts a little smaller
	"legend.fontsize": 8,
	"xtick.labelsize": 8,
	"ytick.labelsize": 8
}

panel_ratio = 0.5
cb_ratio = 0.015

plt.rcParams.update(tex_fonts)
