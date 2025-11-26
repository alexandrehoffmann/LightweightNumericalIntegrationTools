#!/bin/python3

import gnuplot_utils
import matplotlib.pyplot as plt
import math
import numpy as np

import my_plot

path = "../build/demo/"

eps_pow = list(range(4, 16, 2))
sp = (2, len(eps_pow) // 2)

fig, axes = plt.subplots(sp[0], sp[1], sharex=True, sharey=True, figsize=my_plot.get_figure_size(345, subplots=sp), dpi=600, layout='constrained')

lRel_err  = None
lTreshold = None

for ax_nd_id, ax in np.ndenumerate(axes):
	
	ax_id = np.ravel_multi_index(ax_nd_id, sp)
	
	eps = pow(10, -eps_pow[ax_id])
	est_err_mom, est_err = gnuplot_utils.read_2d_curve(path + "addaptive_cch_eps_m" + str(eps_pow[ax_id]) + ".log", 4)
	err_mom,     err     = gnuplot_utils.read_2d_curve(path + "addaptive_cch_eps_m" + str(eps_pow[ax_id]) + ".log", 6)
	
	lRel_err,     = ax.plot(err_mom,     err)
	lRel_est_err, = ax.plot(est_err_mom, est_err)
	
	ax.grid()
	ax.set_yscale('log')
	my_plot.add_label_int(fig, ax, ax_id)

fig.supxlabel("moment order")
fig.legend(handles=[lRel_est_err, lRel_err], labels=["Estimated error", "Real error"], loc='upper center', bbox_to_anchor=(0.5, 1.125), ncol=2)

# ~ plt.savefig('ClensawCurtisErrorEstimate.svg', bbox_inches='tight')
plt.savefig('ClensawCurtisHybridErrorEstimate.pdf', bbox_inches='tight')
plt.clf()
plt.cla()
