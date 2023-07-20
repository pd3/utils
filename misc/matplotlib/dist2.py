#!/usr/bin/python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

import numpy

def smooth(x,window_len=11,window='hanning'):
    if x.ndim != 1: raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len: raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3: return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']: raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s = numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': # moving average
        w = numpy.ones(window_len,'d')
    else:
        w = eval('numpy.'+window+'(window_len)')
    y = numpy.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)]


def bignum(num):
    s = str(num); out = ''; slen = len(s)
    for i in range(slen):
        out += s[i]
        if i+1<slen and (slen-i-1)%3==0: out += ','
    return out

labels  = []
files   = []
files2  = []
# LABELS
# FILES
# FILES2


xdat1 = []
ydat1 = []
xdat2 = []
ydat2 = []
cnts  = []
for i in range(len(files)):
    xdat1.append([])
    ydat1.append([])
    xdat2.append([])
    ydat2.append([])
    fname = files[i]
    with open(fname, 'rb') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if row[0][0] != '#': tmp.append(row)
            else: cnts.append(row[0][2:])
        xdat1[i] = [x[0] for x in tmp]
        ydat1[i] = [x[1] for x in tmp]

    fname = files2[i]
    with open(fname, 'rb') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if row[0][0] != '#': tmp.append(row)
        xdat2[i] = [x[0] for x in tmp]
        ydat2[i] = [float(x[1]) for x in tmp]
        m = max(ydat2[i])
        ydat2[i] = [x/m for x in ydat2[i]]

    
fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(7,7))
for i in range(len(labels)): 
    ax1.plot(xdat1[i],ydat1[i],label='%s (%s)'%(labels[i],bignum(int(cnts[i]))))
for i in range(len(labels)): 
    ax2.plot(xdat2[i],ydat2[i],label='%s (%s)'%(labels[i],bignum(int(cnts[i]))))

ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')
ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')
ax2.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')
ax2.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')

ylabel = None
# yl: ylabel
if ylabel!=None: ax1.set_ylabel(ylabel)
if ylabel!=None: ax2.set_ylabel(ylabel)

xlabel = None
# xl: xlabel
if xlabel!=None: ax2.set_xlabel(xlabel)

xscale = None
# xs: xscale
if xscale!=None: ax1.set_xscale(xscale)     # log, slog
if xscale!=None: ax2.set_xscale(xscale)     # log, slog

title = None
# mt: title
if title!=None: ax1.set_title(title)

plt.legend(numpoints=1,markerscale=2,loc='best',prop={'size':10},frameon=False)
#plt.subplots_adjust(bottom=0.2,left=0.1,right=0.95)

# SAVE
plt.close()

