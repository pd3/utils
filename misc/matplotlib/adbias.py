#!/usr/bin/python

# CMDLINE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv
import numpy
from matplotlib.collections import LineCollection
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

files  = []
# FILES

cols = [ '#337ab7', '#f0ad4e', '#5cb85c', '#5bc0de', '#d9534f', 'grey', 'black' ]

dat = []
for i in range(len(files)):
    fname = files[i]
    with open(fname, 'rb') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if row[0] != 'FT': continue
            af1 = float(row[6])/(float(row[6])+float(row[5]))
            af2 = float(row[8])/(float(row[8])+float(row[7]))
            dat.append([(0,af2),(1,af1)])

wh = (7,5)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)

lt = []
# lt: ['lt']

line_segments = LineCollection(dat,color=cols[0])
ax1.add_collection(line_segments)

#ax1.axes.get_xaxis().set_visible(False)
ax1.set_xticks([0,1])
ax1.set_xticklabels(['fibro','IPS'])
ax1.set_ylim(0,1)

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax1.set_ylabel(ylabel)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel)

xscale = None
# xs: 'xscale'
if xscale!=None: ax1.set_xscale(xscale)     # log, slog

title = None
# title: 'title'
if title!=None: ax1.set_title(title)

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE
plt.close()

