#!/usr/bin/python

# CMDLINE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import csv,sys
import numpy
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

def smooth_data(x,window_len=11,window='hanning'):
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
    odd = 0
    if window_len%2: odd = 1
    y = y[(window_len/2-1):-(window_len/2+odd)]
    return y


def bignum(num):
    s = str(num); out = ''; slen = len(s)
    for i in range(slen):
        out += s[i]
        if i+1<slen and (slen-i-1)%3==0: out += ','
    return out

files  = []
# FILES

if len(files)!=1: sys.exit('Expected one file only for this graph, sorry')
fname = files[0]

(xmin,xmax) = (None,None)
xlim = None
# xr: xlim
if xlim!=None:
    (xmin,xmax) = xlim

xdat = None
ydat = None
cdat = []
with open(fname, 'rb') as f:
    reader = csv.reader(f, 'tab')
    tmp = []
    for row in reader:
        if row[0][0] == '#': continue
        if xmin!=None and xmin>float(row[0]): continue
        if xmax!=None and xmax<float(row[0]): continue
        tmp.append(row)
    xdat = [x[0] for x in tmp]
    ydat = [x[1] for x in tmp]

norm = None     #   +norm 1
# norm: norm
if norm!=None:
    max = float(ydat[0])
    for y in ydat:
        if float(max)<float(y): max = float(y)
    for j in range(len(ydat)):
        ydat[j] = float(ydat[j])/max

sum  = 0
norm = 0
for y in ydat:
    norm += float(y)
for j in range(len(ydat)):
    sum += float(ydat[j])
    cdat.append(sum/norm)

smooth = None
# smooth: smooth
if smooth!=None:
    ydat = smooth_data(numpy.array(ydat,dtype=numpy.float64),smooth,'hanning')


wh = (7,5)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)
ax2 = ax1.twinx()

plt_args1 = {}       # for example: +pa1 "mec='grey',mfc='grey'"
plt_args2 = {}       # for example: +pa2 "mec='grey',mfc='grey'"
# pa1: {plt_args1}
# pa2: {plt_args2}

label1 = 'Cumulative fraction'
label2 = 'Count'
# lb1: 'label1'
# lb2: 'label2'

line1 = '.-'
line2 = '.-'
# lt1: 'line1'
# lt2: 'line2'

col1 = 'black'    # +col1 #4CAE4C (green), #EEA236 (orange), #D43F3A (red), #2E6DA4 (blue)
col2 = '#D43F3A'
# col1: 'col1'
# col2: 'col2'

p1 = ax1.plot(xdat,cdat,line2,label=label1,color=col1,**plt_args1)
p2 = ax2.plot(xdat,ydat,line1,label=label2,color=col2,**plt_args2)
plots = p1 + p2

xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
if xsci!=None: ax1.ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
if ysci!=None: ax2.ticklabel_format(style='sci', scilimits=ysci, axis='y')

grid = None
# gr: grid
if grid!=None: ax1.grid(color='gray', linewidth=1)

if xlim!=None: ax1.set_xlim(xlim)

ylim = None
# yr: 'ylim'
if ylim!=None:                          # for example: +yr 0,1.1%
    (ymin,ymax) = ylim.split(',')
    if ymin[-1] == '%': 
        ymin = float(ymin[0:-1])*min(min(ydat))
    if ymax[-1] == '%': 
        ymax = float(ymax[0:-1])*float(max(max(ydat)))
    ax1.set_ylim(float(ymin),float(ymax))

ax1.set_ylabel(label1, color=col1)
ax2.set_ylabel(label2, color=col2)

for tl in ax1.get_yticklabels(): tl.set_color(col1)
for tl in ax2.get_yticklabels(): tl.set_color(col2)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel)

xscale = None
# xs: 'xscale'
if xscale!=None: ax1.set_xscale(xscale)     # log, slog

yscale = None
# ys: 'yscale'
if yscale!=None: ax2.set_yscale(yscale)     # log, slog

title = None
# title: 'title'
if title!=None: ax1.set_title(title)

labels = [l.get_label() for l in plots]
plt.legend(plots,labels,numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE
plt.close()

