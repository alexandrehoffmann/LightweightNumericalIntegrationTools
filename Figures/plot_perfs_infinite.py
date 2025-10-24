#!/bin/python3

import gbenchmark_utils
import matplotlib.pyplot as plt
import math
import numpy as np

import my_plot

path = "../build/benchmark/"

row_name = "Clensaw Curtis Hybrid"

fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=my_plot.get_figure_size(345, subplots=(1, 1)), dpi=600, layout='constrained')

mom1, cpu_time1 = gbenchmark_utils.read_2d_curve(path + "benchmark.csv", "benchmark_" + row_name.replace(" ", "") + "AdaptiveQuadrature", "cpu_time")
mom2, cpu_time2 = gbenchmark_utils.read_2d_curve(path + "benchmark.csv", "benchmark_Remap" + row_name.replace(" ", "") + "AdaptiveQuadrature", "cpu_time")

ax.plot(mom1, cpu_time1, label="without coordinate remap")
ax.plot(mom2, cpu_time2, label="with coordinate remap")

# ~ print(mom1)
# ~ print(cpu_time1)

ax.set_title(row_name)
ax.set_xlabel("moment order")
ax.set_ylabel("cpu time (s)")
ax.grid()
ax.legend()

fig.tight_layout()
# ~ plt.savefig('perfs_infinite.svg', bbox_inches='tight')
plt.savefig('perfs_infinite.pdf', bbox_inches='tight')
plt.clf()
plt.cla()
