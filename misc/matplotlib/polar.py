#!/usr/bin/python
#
#
#   see also ~/wtxt/logs/sandbox/other/polar-mri.py
#
#

# CMDLINE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv,sys
import numpy as np
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

files  = []
# FILES

# mpl.rcParams['axes.color_cycle'] = [ '#337ab7', '#f0ad4e', '#5cb85c', '#5bc0de', '#d9534f', 'grey', 'black' ]

lines = open(files[0], 'rb').readlines()

if len(lines) != 2: 
    print "Error: expected two lines on input, got: ",len(lines)
    sys.exit(1)

cols  = lines[0].rstrip('\n').split('\t')
yvals = lines[1].rstrip('\n').split('\t')

if len(cols) != len(yvals): 
    print "Error: different number of fields in the two lines: ",len(cols)," vs ",len(yvals)
    sys.exit(1)

cols  = cols[4:]
yvals = yvals[4:]
xvals = []
for i in range(len(yvals)):
    xvals.append(i*2*np.pi/len(yvals))
    yvals[i] = float(yvals[i].rstrip('%')) + 100.

width = 2*np.pi / len(yvals)

ax = plt.subplot(111, projection='polar')
bars = ax.bar(xvals, yvals, width=width, bottom=0.0)

for x,y,bar in zip(xvals, yvals, bars):
    if abs(y-100) > 5:
        bar.set_facecolor('red')
        bar.set_edgecolor('red')
    else:
        bar.set_facecolor('grey')
        bar.set_edgecolor('grey')
    bar.set_alpha(0.5)


# dpi: dpi

# SAVE
plt.close()

