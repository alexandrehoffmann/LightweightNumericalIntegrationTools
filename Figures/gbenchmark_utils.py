## module gnuplot_utils## module gbenchmark_utils

import csv

def read_2d_curve(fname, row_label, col_label="cpu_time"):
	time_mulitplier = {
		"ns": 1e-9,
		"us": 1e-6,
		"ms": 1e-3
	}
	
	x = []
	fx = []
	with open(fname) as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',', quotechar='"')
		first_tokens = next(spamreader)
		i = first_tokens.index(col_label)
		for tokens in spamreader:
			if row_label in tokens[0]:
				x.append(float(tokens[0].split('/')[1]))
				fx.append(float(tokens[i])*time_mulitplier[tokens[i+1]])
	return x, fx

