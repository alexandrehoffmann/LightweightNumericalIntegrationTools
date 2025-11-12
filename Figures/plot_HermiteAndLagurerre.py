#!/bin/python3

import gnuplot_utils
import matplotlib.pyplot as plt
import math
import numpy as np

import my_plot

path = "../build/demo/"

fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=my_plot.get_figure_size(345, subplots=(1, 1)), dpi=600, layout='constrained')

mom, herite_rel_err = gnuplot_utils.read_2d_curve(path + "hermite_and_laguerre.log", 6)
mom, lagurre_rel_err = gnuplot_utils.read_2d_curve(path + "hermite_and_laguerre.log", 7)

ax.plot(mom[0:40], herite_rel_err[0:40], label="Guass-Hermite")
ax.plot(mom[0:40], lagurre_rel_err[0:40], label="Guass-Laguerre")

ax.set_xlabel("moment order")
ax.set_ylabel("relative error")
ax.grid()
ax.legend()

fig.tight_layout()
plt.savefig('HermiteAndLagurerre.svg', bbox_inches='tight')
# ~ plt.savefig('HermiteAndLagurerre.pdf', bbox_inches='tight')
plt.clf()
plt.cla()
