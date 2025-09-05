#!/bin/python3

import gnuplot_utils
import matplotlib.pyplot as plt
import math
import numpy as np

import my_plot

path = "../build/demo/"

eps_pow = list(range(4, 16, 2))
sp = (2, len(eps_pow) // 2)

fig, axes = plt.subplots(sp[0], sp[1], sharex=True, sharey=True, figsize=my_plot.get_figure_size(2*345, subplots=sp), dpi=600, layout='constrained')

lRel_err  = None
lTreshold = None

for ax_nd_id, ax in np.ndenumerate(axes):
	
	ax_id = np.ravel_multi_index(ax_nd_id, sp)
	
	eps = pow(10, -eps_pow[ax_id])
	mom, rel_err = gnuplot_utils.read_2d_curve(path + "addaptive_cch_eps_m" + str(eps_pow[ax_id]) + ".log", 6)
	
	lRel_err  = ax.plot(mom, rel_err)
	lTreshold = ax.plot(mom, eps*np.ones(len(mom)))
	
	ax.grid()
	ax.set_yscale('log')
	my_plot.add_label_int(fig, ax, ax_id)

fig.supxlabel("moment order")

bb = (fig.subplotpars.left, fig.subplotpars.top, fig.subplotpars.right-fig.subplotpars.left,.065)

leg = axes[0, 0].legend(
	[lRel_err, lTreshold], 
	labels=["Clensaw Curtis relative error", "treshold"], 
	ncol=2, 
	bbox_to_anchor=bb, 
	mode="expand", 
	borderaxespad=0., 
	bbox_transform=fig.transFigure)

fig.tight_layout()
fig.subplots_adjust(hspace=0.04, wspace=0.04)

plt.savefig('ClensawCurtisHybridError.pdf', bbox_inches='tight')
plt.clf()
plt.cla()
