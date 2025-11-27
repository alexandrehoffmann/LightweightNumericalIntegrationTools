#!/bin/python3

import gbenchmark_utils
import matplotlib.pyplot as plt
import math
import numpy as np

import my_plot

path = "../build/benchmark/"

row_names = ["Gauss Legendre", "Clensaw Curtis", "Clensaw Curtis Hybrid", "GL CC"]

fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=my_plot.get_figure_size(345, subplots=(1, 1)), dpi=600, layout='constrained')

for row_name in row_names:
	mom, cpu_time = gbenchmark_utils.read_2d_curve(path + "benchmark.csv", "benchmark_Remap" + row_name.replace(" ", "") + "AdaptiveQuadrature", "cpu_time")
	ax.plot(mom, cpu_time, label=row_name)

ax.set_xlabel("moment order")
ax.set_ylabel("cpu time (s)")
ax.grid()
ax.legend()

fig.tight_layout()
plt.savefig('perfs_remap.svg', bbox_inches='tight')
# ~ plt.savefig('perfs_remap.pdf', bbox_inches='tight')
plt.clf()
plt.cla()
