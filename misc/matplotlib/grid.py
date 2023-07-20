#!/usr/bin/python

# CMDLINE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv
import numpy
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

nrc = (1,1)     # nrows,ncols
# nrc: nrc

style = None        # xkcd, ggplot, ...
# sty: 'style'
if style=='xkcd':
    plt.xkcd()
elif style!=None:
    plt.style.use(style)

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

(xmin,xmax) = (None,None)
xlim = None
# xr: xlim
if xlim!=None:
    (xmin,xmax) = xlim

mpl.rcParams['axes.color_cycle'] = [ '#337ab7', '#f0ad4e', '#5cb85c', '#5bc0de', '#d9534f', 'grey', 'black' ]

xdat = []
ydat = []
for i in range(len(files)):
    xdat.append([])
    ydat.append([])
    fname = files[i]
    with open(fname, 'rb') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if row[0][0] == '#': continue
            if xmin!=None and xmin>float(row[0]): continue
            if xmax!=None and xmax<float(row[0]): continue
            tmp.append(row)
        xdat[i] = [x[0] for x in tmp]
        ydat[i] = [x[1] for x in tmp]

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

fig, ax = plt.subplots(nrc[0], nrc[1], figsize=wh)
if type(ax).__name__!='ndarray': ax = numpy.array([ax],ndmin=2)
elif type(ax[0]).__name__!='ndarray': ax = numpy.array([ax],ndmin=2)
ax1 = ax[0][0]

plt_args = {}       # for example: +pa "mec='grey',mfc='grey'"
# pa: {plt_args}

smooth = None
# smooth: smooth
lb = []
# lb: ['lb']
lt = []
# lt: ['lt']
xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
grid = None
# gr: grid
ylim = None
# yr: 'ylim'
ylabel = None
# yl: 'ylabel'
xlabel = None
# xl: 'xlabel'
xscale = None
# xs: 'xscale'
yscale = None
# ys: 'yscale'

i = 0
for irow in range(len(ax)):
    for icol in range(len(ax[irow])):
        line = '-'
        if i < len(lt): line = lt[i]
        label = ''
        if i <len(lb): label = lb[i]
        if smooth!=None:
            ydat[i] = smooth_data(numpy.array(ydat[i],dtype=numpy.float64),smooth,'hanning')
        ax[irow][icol].plot(xdat[i],ydat[i],line,label=label,**plt_args)
        i = i + 1

        if xsci!=None: ax[irow][icol].ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
        if ysci!=None: ax[irow][icol].ticklabel_format(style='sci', scilimits=ysci, axis='y')
        if grid!=None: ax[irow][icol].grid(color='gray', linewidth=1)
        if xlim!=None: ax[irow][icol].set_xlim(xlim)
        if ylim!=None:                          # for example: +yr 0,1.1%
            (ymin,ymax) = ylim.split(',')
            if ymin[-1] == '%': 
                ymin = float(ymin[0:-1])*min(min(ydat))
            if ymax[-1] == '%': 
                ymax = float(ymax[0:-1])*float(max(max(ydat)))
            ax[irow][icol].set_ylim(float(ymin),float(ymax))
        
        if ylabel!=None: ax[irow][icol].set_ylabel(ylabel)
        if xlabel!=None: ax[irow][icol].set_xlabel(xlabel)
        if xscale!=None: ax[irow][icol].set_xscale(xscale)     # log, symlog
        if yscale!=None: ax[irow][icol].set_yscale(yscale)     # log, symlog
        args = {}
        if style=='ggplot':
            ax[irow][icol].spines['top'].set_visible(False)
            ax[irow][icol].spines['right'].set_visible(False)
            ax[irow][icol].get_xaxis().tick_bottom()
            ax[irow][icol].get_yaxis().tick_left()
            ax[irow][icol].spines['bottom'].set_color('grey')
            ax[irow][icol].spines['left'].set_color('grey')
            mpl.rcParams['text.color'] = '555555'
            args = {'color':'#555555'}

title = None
# title: 'title'
if title!=None: ax1.set_title(title,**args)

if len(lb)>0: plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE
plt.close()

