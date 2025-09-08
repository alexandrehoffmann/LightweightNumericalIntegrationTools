#!/bin/python3

import gbenchmark_utils
import matplotlib.pyplot as plt
import math
import numpy as np

import my_plot

path = "../build/benchmark/"

row_names = ["Gauss Legendre", "Clensaw Curtis", "Clensaw Curtis Hybrid"]

fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=my_plot.get_figure_size(2*345, subplots=(1, 1)), dpi=600, layout='constrained')

for row_name in row_names:
	mom, cpu_time = gbenchmark_utils.read_2d_curve(path + "benchmark.csv", "benchmark_" + row_name.replace(" ", "") + "AdaptiveQuadrature", "cpu_time")
	ax.plot(mom, cpu_time, label=row_name)

ax.set_xlabel("moment order")
ax.set_ylabel("cpu time (s)")
ax.grid()
ax.legend()

fig.tight_layout()
plt.savefig('perfs.svg', bbox_inches='tight')
plt.clf()
plt.cla()
