## module gnuplot_utils

def read_2d_curve(fname, i=1):
	x = []
	fx = []
	with open(fname) as f:
		for _line in f.readlines():
			line = _line.strip()
			if line[0] != '#':
				tokens =  [float(token.replace("false", "0").replace("true", "1")) for token in filter(None, line.split(" "))]
				x.append(tokens[0])
				fx.append( tokens[i] )
	return x, fx
