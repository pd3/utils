#!/usr/bin/python

# CMDLINE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv
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

style = None        # mine, xkcd, ggplot, ...
# sty: 'style'
if style=='xkcd':
    plt.xkcd()
elif style!=None and style!='mine':
    plt.style.use(style)

files  = []
# FILES

mpl.rcParams['axes.color_cycle'] = [ '#337ab7', '#f0ad4e', '#5cb85c', '#5bc0de', '#d9534f', 'grey', 'black' ]

ydat = []
cnts = []
for i in range(len(files)):
    ydat.append([])
    fname = files[i]
    with open(fname, 'rb') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if row[0][0] != '#': tmp.append(row)
            else: cnts.append(row[0][2:])
        ydat[i] = [x[0] for x in tmp]

norm = None
# norm: norm
if norm!=None:
    for i in range(len(files)):
        max = float(ydat[i][0])
        for y in ydat[i]:
            if float(max)<float(y): max = float(y)
        for j in range(len(ydat[i])):
            ydat[i][j] = float(ydat[i][j])/max


wh = (7,5)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)

plt_args = {}       # for example: +pa "mec='grey',mfc='grey',alpha=0.7"
# pa: {plt_args}

smooth = None
# smooth: smooth
lb = []
# lb: ['lb']
lt = []
# lt: ['lt']
cl = None
# cl: ['cl']
for i in range(len(files)):
    line = '-'
    if i < len(lt): line = lt[i]
    label = ''
    if i <len(lb): label = lb[i]
    if cl!=None: plt_args['color'] = cl[i]
    if smooth!=None:
        ydat[i] = smooth_data(numpy.array(ydat[i],dtype=numpy.float64),smooth,'hanning')
    ax1.plot(ydat[i],line,label=label,**plt_args)

xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
if xsci!=None: ax1.ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
if ysci!=None: ax1.ticklabel_format(style='sci', scilimits=ysci, axis='y')

xticks = None       # +xticks -
# xticks: 'xticks'
if xticks=='-': ax1.set_xticks([])

yticks = None       # +yticks -
# yticks: 'yticks'
if yticks=='-': ax1.set_yticks([])

grid = None
# gr: grid
if grid!=None: ax1.grid(color='gray', linewidth=1)

xlim = None
# xr: 'xlim'
if xlim!=None:                          # for example: +xr -10%,110%
    (xmin,xmax) = xlim.split(',')
    if xmin[-1] == '%': 
        xmin = 0.01*float(xmin[0:-1])*len(ydat[0])
    if xmax[-1] == '%': 
        xmax = 0.01*float(xmax[0:-1])*len(ydat[0])
    ax1.set_xlim(float(xmin),float(xmax))

ylim = None
# yr: ylim
if ylim!=None: ax1.set_ylim(ylim)

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax1.set_ylabel(ylabel)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel)

xscale = None
# xs: 'xscale'
if xscale!=None: ax1.set_xscale(xscale)     # log, slog

args = {}
if style=='ggplot':
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('grey')
    ax1.spines['left'].set_color('grey')
    mpl.rcParams['text.color'] = '555555'
    args = {'color':'#555555'}
if style=='mine':
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('grey')
    ax1.spines['left'].set_color('grey')
    mpl.rcParams['text.color'] = '555555'
    args = {'color':'#555555'}
    ax1.patch.set_visible(False)

title = None
# title: 'title'
if title!=None: ax1.set_title(title)

if len(lb)>0: plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE
plt.close()

