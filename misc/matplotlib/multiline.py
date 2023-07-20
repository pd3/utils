#!/usr/bin/python
#
#   x1  y11 y12 y13 y14 .. y1n
#   ..  ..
#   xm  ym1 ym2 ym3 ym4 .. ymn
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import itertools
import csv
import math
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

xdat   = []
dat    = None
files  = []
labels = []
# FILES
# LABELS

with open(files[0], 'r') as f:
    reader = csv.reader(f, 'tab')
    for row in reader:
        if row[0][0] == '#': continue
        if dat==None:
            dat = []
            for i in range(1,len(row)): dat.append([])
        for i in range(1,len(row)):
            dat[i-1].append(row[i])
        xdat.append(row[0])

wh = (7,5)
# wh: wh
fig, ax = plt.subplots(1, 1, figsize=wh)

ylabel = None
# yl: 'ylabel'

xlabel = None
# xl: 'xlabel'

title = None
# title: 'title'
if title!=None: ax.set_title(title)

labelsize = None
# ls: labelsize

for i in range(len(dat)):
    ax.plot(xdat,dat[i],'-',color='gray')

xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
if xsci!=None: ax.ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
if ysci!=None: ax.ticklabel_format(style='sci', scilimits=ysci, axis='y')

#plt.subplots_adjust(bottom=0.05,left=0.05,right=0.95,top=0.95,hspace=0.12,wspace=0.12)

# SAVE
plt.close()

