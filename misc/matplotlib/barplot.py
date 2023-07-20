#!/usr/bin/env python3
#
#   cat cnt.txt     | mplot barplot -o barplot.png +type cnt
#   cat x-cnt.txt   | mplot barplot -o barplot.png +type x-cnt
#   cat lb-cnt.txt  | mplot barplot -o barplot.png +type lbl-cnt
#   cat hist.txt    | mplot barplot -o barplot.png +type hist 
#   cat ylb-cnt.txt | mplot barplot -o barplot.png +type ylbl-cnt
#   cat lb-cnt-color.txt | mplot barplot -o barplot.png +type lbl-cnt-col
#
#   Options
#       draw a line
#           +line '42,0,42,1e9;col="red";lw=2;lt=":"'

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv,re
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

style = None        # xkcd, ggplot, ...
# sty: 'style'
if style=='xkcd':
    plt.xkcd()
elif style!=None and style!='mine':
    plt.style.use(style)

def smooth_data(x,window_len=11,window='hanning'):
    if x.ndim != 1: raise Exception("smooth only accepts 1 dimension arrays.")
    if x.size < window_len: raise Exception("Input vector needs to be bigger than window size.")
    if window_len<3: return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']: raise Exception("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s = numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': # moving average
        w = numpy.ones(window_len,'d')
    else:
        w = eval('numpy.'+window+'(window_len)')
    y = numpy.convolve(w/w.sum(),s,mode='valid')
    odd = 0
    if window_len%2: odd = 1
    y = y[(int(window_len/2)-1):-(int(window_len/2)+odd)]
    return y

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

type = 'cnt'
# type: 'type'

direction = 'vertical'
if type=='ylbl-cnt':
    direction = 'horizontal'

lb = []
# lb: ['lb']

xdat = []
ydat = []
lbls = []
cols = []
for i in range(len(files)):
    xdat.append([])
    ydat.append([])
    lbls.append([])
    cols.append([])
    fname = files[i]
    with open(fname, 'r') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if row[0][0] == '#': continue
            if type=='cnt':
                ydat[i].append(float(row[0]))
                xdat[i].append(len(ydat[i]))
            elif type=='lbl-cnt':
                row[0] = re.sub(r'\\n',"\n",row[0])
                lbls[i].append(row[0])
                ydat[i].append(float(row[1]))
                xdat[i].append(len(ydat[i]))
            elif type=='lbl-cnt-col':
                row[0] = re.sub(r'\\n',"\n",row[0])
                lbls[i].append(row[0])
                ydat[i].append(float(row[1]))
                xdat[i].append(len(ydat[i]))
                cols[i].append(row[2])
            elif type=='ylbl-cnt':
                row[0] = re.sub(r'\\n',"\n",row[0])
                lbls[i].append(row[0])
                ydat[i].append(float(row[1]))
                xdat[i].append(len(ydat[i]))
            elif type=='x-cnt':
                xdat[i].append(float(row[0]))
                ydat[i].append(float(row[1]))
            elif type=='hist':
                ydat[i].append(float(row[0]))
            else:
                print >> sys.stderr, "Not supported: +type ",type
                sys.exit()

smooth = None       # +smooth 3
# smooth: smooth
if smooth!=None: import numpy

wh = (7,5)
# wh: wh

wd = 0.7        # fraction of dx distance
# wd: wd

min_dx = None
for dat in xdat:
    for i in range(len(dat)-1):
        if min_dx==None or min_dx > abs(dat[i+1]-dat[i]): min_dx = abs(dat[i+1]-dat[i])
if min_dx==None: min_dx = 1
wd = min_dx*wd

fcolor = None
# fc: 'fcolor'

ecolor = None
# ec: 'ecolor'

colors = [ '#337ab7', '#f0ad4e', '#5cb85c', '#5bc0de', '#d9534f', 'grey', 'black' ]

plt_args = {}       # for example: +pa "mec='grey',mfc='grey'"
# pa: {plt_args}

if direction=='vertical':
    lbl_args = {'rotation':35,'ha':'right','multialignment':'center'} # for example: "rotation=35,ha='right',ma='center',fontsize=9"
else:
    lbl_args = {'ha':'right','multialignment':'center'} # for example: "rotation=35,ha='right',ma='center',fontsize=9"
# la: {lbl_args}

fig, ax1 = plt.subplots(1, 1, figsize=wh)
for i in range(len(labels)): 
    plt_args['color'] = colors[i%len(colors)]
    if fcolor!=None: plt_args['color'] = fcolor
    if ecolor==None: plt_args['edgecolor'] = plt_args['color']
    if len(cols[i]):
        plt_args['color'] = cols[i]
        if ecolor==None: plt_args['edgecolor'] = cols[i]
    if i < len(lb): plt_args['label'] = lb[i]
    if smooth!=None:
        ydat[i] = smooth_data(numpy.array(ydat[i],dtype=numpy.float64),smooth,'hanning')
    if type=='hist':
        ax1.hist(ydat[i],wd,**plt_args)
    elif direction=='horizontal':
        ax1.barh([x for x in xdat[i]],ydat[i],wd,**plt_args)
        if len(lbls[i]):
            ax1.set_yticks([x for x in xdat[i]])
            ax1.set_yticklabels(lbls[i],**lbl_args)
    else:
        ax1.bar([x for x in xdat[i]],ydat[i],wd,**plt_args)
        if len(lbls[i]):
            ax1.set_xticks([x for x in xdat[i]])
            ax1.set_xticklabels(lbls[i],**lbl_args)


line = None         # +line '42,0,42,1e9;col="red";lw=2;ls=":"'
# line: 'line'
if line!=None:
    tmp = line.split(';')
    lin = tmp[0].split(',')
    tmp_args = {}
    for x in tmp[1:]:
        a,b = x.split('=')
        tmp_args[a] = b
    ax1.plot([float(lin[0]),float(lin[2])],[float(lin[1]),float(lin[3])],**tmp_args)

hdt = None       # hide ticks: +hdt 1
# hdt: hdt
if hdt:
    if direction=='horizontal':
        ax1.yaxis.set_tick_params(length=0)

if style=='ggplot':
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('grey')
    ax1.spines['left'].set_color('grey')
    mpl.rcParams['text.color'] = '555555'
    args = {'color':'#555555'}
if style=='mine' or style=='xkcd':
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('grey')
    ax1.spines['left'].set_color('grey')
    mpl.rcParams['text.color'] = '555555'
    #args = {'color':'#555555'}
    ax1.patch.set_visible(False)

xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
if xsci!=None: ax1.ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
if ysci!=None: ax1.ticklabel_format(style='sci', scilimits=ysci, axis='y')

xlim = None
# xr: 'xlim'
if xlim!=None:
    x = xlim.split(',')
    xlim = {}
    if x[0]!='': xlim['left'] = float(x[0])
    if x[1]!='': xlim['right'] = float(x[1])
    ax1.set_xlim(**xlim)

ylim = None
# yr: 'ylim'
if ylim!=None:
    y = ylim.split(',')
    ylim = {}
    if y[0]!='': ylim['bottom'] = float(y[0])
    if y[1]!='': ylim['top'] = float(y[1])
    ax1.set_ylim(**ylim)

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax1.set_ylabel(ylabel,labelpad=10)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel,labelpad=10)

xticks = 1       # +xt 0
# xt: xticks
if __builtins__.type(xticks).__name__ == 'tuple': ax1.set_xticks(xticks)
elif xticks==0: ax1.set_xticks([])

xscale = None
# xs: 'xscale'
if xscale!=None: ax1.set_xscale(xscale)     # log, slog

yscale = None
# ys: 'yscale'
if yscale!=None: ax1.set_yscale(yscale)     # log, slog

title = None
# title: 'title'
if title!=None: ax1.set_title(title,y=1.08)


if len(lb)>0: plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE
plt.close()

