#!/usr/bin/python

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

labels = []
files  = []
# LABELS
# FILES


xdat = []
ydat = []
cnts = []
for i in range(len(files)):
    xdat.append([])
    ydat.append([])
    fname = files[i]
    with open(fname, 'rb') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if row[0][0] != '#': tmp.append(row)
            else: cnts.append(row[0][2:])
        xdat[i] = [x[0] for x in tmp]
        ydat[i] = [x[1] for x in tmp]

fig, ax1 = plt.subplots(1, 1, figsize=(7,5))
for i in range(len(labels)): 
    ax1.plot(xdat[i],ydat[i],label='%s (%s)'%(labels[i],bignum(int(cnts[i]))))

ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')
ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')

ylabel = None
# yl: ylabel
if ylabel!=None: ax1.set_ylabel(ylabel)

xlabel = None
# xl: xlabel
if xlabel!=None: ax1.set_xlabel(xlabel)

xscale = None
# xs: xscale
if xscale!=None: ax1.set_xscale(xscale)     # log, slog

plt.legend(numpoints=1,markerscale=2,loc='best',prop={'size':10},frameon=False)
#plt.subplots_adjust(bottom=0.2,left=0.1,right=0.95)

# SAVE
plt.close()

