#!/usr/bin/python
#
# Creates multiple bar plots on one page
#
# cat plots/pophist.txt | head -3 | cut -f5- | mplot bartable -o plots/pophist.png +lb CN3,CN4,CN5,CN6,CN7
#
# - arbitrary number of rows and columns
# - rows: single dataset (e.g. a sample)
# - columns: categories
#
# For example:
#   37163   24809   3   0   201 256
#   32405   27151   0   143 593 719
#   33441   24743   10  10  265 466
#
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
            xdat.append(len(ydat[0])-0.5)
            for icol in range(len(ydat)):
                ydat[icol].append(float(row[icol]))

labels = []
# lb: ('labels')

if len(labels) and len(labels)!=len(ydat):
    print "Wrong number of labels: ",len(labels),len(ydat)

wh = (7,5)
# wh: wh

ysci = None
# ysci: (ysci)  # +ysci -2,2

ecolor = None
# ecolor: ec

colors = [ '#337ab7', '#f0ad4e', '#5cb85c', '#5bc0de', '#d9534f', 'grey', 'black' ]

fig, ax = plt.subplots(len(ydat), 1, figsize=wh, sharex=True)
for i in range(len(ydat)): 
    args = {}
    args['color'] = colors[i%len(colors)]
    if ecolor==None: args['edgecolor'] = args['color']
    ax[i].bar(xdat, ydat[i], 1, **args)
    if ysci!=None: ax[i].ticklabel_format(style='sci', scilimits=ysci, axis='y')
    if len(labels): ax[i].set_ylabel(labels[i])

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax[-1].set_xlabel(xlabel)

xlim = None
# xr: xlim
if xlim!=None: ax[-1].set_xlim(xlim)

adjust = {'left':0.08,'bottom':0.08,'right':0.95,'top':0.95}
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# SAVE
plt.close()

