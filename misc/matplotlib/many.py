#!/usr/bin/env python3
#
#   cat dat.txt   | mplot many -o img.png
#
#   dat.txt
#       1st row .. x-values (colors, see below)
#       2+  row .. y-values of curves to plot
#
#   Colors
#       Each line can be colored separately if the first value of the first line is 'color'.
#       Note that color hash strings must be escaped with '\', otherwise the line would be
#       removed by mplot.
#
# CMDLINE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys

style = None        # xkcd, ggplot, ...
# sty: 'style'
if style=='xkcd':
    plt.xkcd()
elif style!=None and style!='mine':
    plt.style.use(style)


files  = []
# FILES
fname = files[0]


color = 'grey'      # +cl red
# cl: 'color'


(xmin,xmax) = (None,None)
xlim = None         # +xr 1,10    +xr 1,    +xr ,10
# xr: 'xlim'
if xlim!=None:
    x = xlim.split(',')
    xlim = {}
    if x[0]!='': xmin = float(x[0]); xlim['left']  = xmin
    if x[1]!='': xmax = float(x[1]); xlim['right'] = xmax

colors = None
xrow = None
dat  = []
reader = open(fname)
for line in reader:
    row = line.rstrip('\n').split('\t')
    if xrow==None:
        if row[0]=='color':
            colors = []
            xrow = row[1:]
        else:
            xrow = row
        for i in range(len(xrow)): xrow[i] = float(xrow[i])
        continue
    if colors!=None:
        if row[0][0]=='\\': row[0] = row[0][1:]
        colors.append(row[0])
        row = row[1:]
    for i in range(len(row)): row[i] = float(row[i])
    dat.append([xrow,row])

yscale = None       # +ys log,symlog
# ys: 'yscale'

xscale = None       # +xs log,symlog
# xs: 'xscale'

wh = (7,5)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)

asp = None  # +asp 1
# asp: asp
if asp != None: ax1.set_aspect(asp, adjustable='box')

plt_args = {}       # for example: +pa "mec='grey',mfc='grey',zorder=100,clip_on=False,alpha=0.5"
if style=='mine': plt_args = {'zorder':100,'clip_on':False}
# pa: {plt_args}

lb = []
# lb: ['lb']
lt = []             # line types:   +lt --
# lt: ['lt']
lc = []             # line colors:  +lc '#337ab7'
# lc: ['lc']
for i in range(len(dat)):
    row = dat[i]
    line = '-'
    if colors==None:
        plt_args['color'] = color
    else:
        plt_args['color'] = colors[i]
    ax1.plot(row[0],row[1],line,**plt_args)

xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
if xsci!=None: ax1.ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
if ysci!=None: ax1.ticklabel_format(style='sci', scilimits=ysci, axis='y')

xt = None       # xticks: +xt 1,2,3
# xt: 'xt'
if xt!=None:
    tmp = xt.split(',')
    ax1.set_xticks([float(x) for x in tmp])

xlab_args = {}      # for example: +xarg "rotation=45,ha='right',ma='left',fontsize=9"
# xarg: {xlab_args}         

xtl = None      # xticklabels: +xtl 'label 1;label 2'
# xtl: 'xtl'
if xtl!=None:
    #ax1.get_xaxis().set_tick_params(direction='out')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xticklabels(xtl.split(';'), **xlab_args)

grid = None     # +gr '#eeeeee,1,--'    # color,linewidth,linestyle
# gr: 'grid'
if grid!=None:
    grid = grid.split(',')
    ax1.grid(color=grid[0], linewidth=float(grid[1]), linestyle=grid[2])

if xlim!=None: ax1.set_xlim(**xlim)

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
if ylabel!=None: ax1.set_ylabel(ylabel,labelpad=10)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel,labelpad=10)

if yscale!=None: ax1.set_yscale(yscale)     # log, symlog
if xscale!=None: ax1.set_xscale(xscale)     # log, symlog

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
    #args = {'color':'#555555'}
    ax1.patch.set_visible(False)

ta_args = {}        #  +ta t=1.0
# ta: {ta_args}
ta_args = {**args,**ta_args}
title = None
# title: 'title'
if title!=None: ax1.set_title(title,**ta_args)

lg_args = {}            # +lga title='Quality Decile'
# lga: {lg_args}
lg_args = dict({'numpoints':1,'markerscale':1,'loc':'best','prop':{'size':10},'frameon':False},**lg_args)

if len(lb)>0: plt.legend(**lg_args)

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

