#!/usr/bin/python
#
# cat plots/pophist.txt | head -3 | cut -f5- | mplot bartable -o plots/pophist.png +lb CN3,CN4,CN5,CN6,CN7
#
# - arbitrary number of rows and columns
# - each row corresponds to a single bar cluster
# - first column is displayed as a first bar in each cluster
#
# For example:
#   37163   24809   3   0   201 256
#   32405   27151   0   143 593 719
#   33441   24743   10  10  265 466
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

def bignum(num):
    s = str(num); out = ''; slen = len(s)
    for i in range(slen):
        out += s[i]
        if i+1<slen and (slen-i-1)%3==0: out += ','
    return out

files  = []
# FILES


xdat = []   # location of the groups
ydat = []
for i in range(len(files)):
    fname = files[i]
    with open(fname, 'rb') as f:
        reader = csv.reader(f, 'tab')
        for row in reader:
            if row[0][0] == '#': continue
            if len(ydat)==0:
                for x in row: 
                    ydat.append([])
            xdat.append(len(ydat[0]))
            for icol in range(len(ydat)):
                ydat[icol].append(float(row[icol]))

labels = []
# lb: ('labels')

if len(labels) and len(labels)!=len(ydat):
    print "Wrong number of labels: ",len(labels),len(ydat)

wh = (7,5)
# wh: wh

wd = 1./(2 + len(ydat))
# wd: wd

ecolor = None
# ec: 'ecolor'

colors = [ '#337ab7', '#f0ad4e', '#5cb85c', '#5bc0de', '#d9534f', 'grey', 'black' ]

fig, ax1 = plt.subplots(1, 1, figsize=wh)
for i in range(len(ydat)): 
    args = {}
    args['color'] = colors[i%len(colors)]
    if len(labels): args['label'] = labels[i]
    if ecolor==None: args['edgecolor'] = args['color']
    ax1.bar([x+i*wd for x in xdat], ydat[i], wd, **args)

#ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')
#ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')


xtick_args = {} # for example: +xta "rotation=35,ha='right'"
# xta: {xtick_args}

xticks = None
# xt: ('xticks')
if xticks!=None:
    ax1.xaxis.set_ticks_position('none')
    ticks = range(len(ydat[0]))
    ax1.set_xticks([x+0.5 for x in ticks])
    ax1.set_xticklabels(xticks,**xtick_args)

ax1.set_xlim(-wd,len(ydat[0])+wd)

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax1.set_ylabel(ylabel)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel)

xscale = None
# xs: xscale
if xscale!=None: ax1.set_xscale(xscale)     # log, slog

title = None
# title: 'title'
if title!=None: ax1.set_title(title)

if len(labels)>0: plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':8},frameon=False,ncol=4)

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: +adj bottom=0.2,left=0.1,right=0.95

# SAVE
plt.close()

