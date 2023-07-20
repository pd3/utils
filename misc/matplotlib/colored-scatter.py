#!/usr/bin/python
#
#   x y value
#

# CMDLINE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv
import numpy
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

style = None        # xkcd, ggplot, ...
# sty: 'style'
if style=='xkcd':
    plt.xkcd()
elif style!=None and style!='mine':
    plt.style.use(style)

files  = []
# FILES

if len(files)!=1: print 'Expected one file only'; sys.exit(1)

(xmin,xmax) = (None,None)
xlim = None
# xr: xlim
if xlim!=None:
    (xmin,xmax) = xlim

mpl.rcParams['axes.color_cycle'] = [ '#337ab7', '#f0ad4e', '#5cb85c', '#5bc0de', '#d9534f', 'grey', 'black' ]

xdat = []
ydat = []
zdat = []
for i in range(len(files)):
    xdat.append([])
    ydat.append([])
    zdat.append([])
    fname = files[i]
    with open(fname, 'r') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if row[0][0] == '#': continue
            if xmin!=None and xmin>float(row[0]): continue
            if xmax!=None and xmax<float(row[0]): continue
            tmp.append(row)
        xdat[i] = [float(x[0]) for x in tmp]
        ydat[i] = [float(x[1]) for x in tmp]
        zdat[i] = [float(x[2]) for x in tmp]

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

asp = None  # +asp 1
# asp: asp
if asp != None: ax1.set_aspect(asp, adjustable='box')

plt_args = {}       # for example: +pa "mec='grey',mfc='grey'"
# pa: {plt_args}

line = None     # eg. +line '0,1,0,1,k--'
# line: 'line'
if line!=None:
    x = line.split(',')
    ax1.plot([x[0],x[1]],[x[2],x[3]],x[4])

cscale = None
# cs: 'cscale'
if cscale!=None: plt_args['norm'] = mpl.colors.LogNorm()


smooth = None
# smooth: smooth
lb = []
# lb: ['lb']
lt = []
# lt: ['lt']
for i in range(len(files)):
    line = '-'
    if i < len(lt): line = lt[i]
    label = ''
    if i <len(lb): label = lb[i]
    sc = ax1.scatter(xdat[i],ydat[i], edgecolors='none', marker='o', c=zdat[i], **plt_args)

xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
if xsci!=None: ax1.ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
if ysci!=None: ax1.ticklabel_format(style='sci', scilimits=ysci, axis='y')

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

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax1.set_ylabel(ylabel)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel)

clabel = None
# cl: 'clabel'
if clabel!=None:
    cb = fig.colorbar(sc, ax=ax1, format='%.2e');
    cb.set_label(clabel)
    for t in cb.ax.get_yticklabels(): t.set_fontsize(9)

xscale = None
# xs: 'xscale'
if xscale!=None: ax1.set_xscale(xscale)     # log, symlog

yscale = None
# ys: 'yscale'
if yscale!=None: ax1.set_yscale(yscale)     # log, symlog

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
if title!=None: ax1.set_title(title,**args)

if len(lb)>0: plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)

rect = None
# rect: rect
if rect!=None:
    import matplotlib.patches as patches
    y0,y1 = ax1.get_ylim()
    width = rect[1]-rect[0]
    xy = [rect[0],y0]
    height = y1-y0 
    rect = patches.Rectangle(xy, width, height, color='#F2DEDE',zorder=-100)
    ax1.add_patch(rect)

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE
plt.close()

