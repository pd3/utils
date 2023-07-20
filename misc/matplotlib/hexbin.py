#!/usr/bin/env python3
#
#   x y
#

# CMDLINE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

style = None        # xkcd, ggplot, ...
# sty: 'style'
if style=='xkcd':
    plt.xkcd()
elif style!=None and style!='mine':
    plt.style.use(style)

files  = []
# FILES

cdat = None
# cdat: cdat
if cdat!=None: cdat = []

xdat = []
ydat = []
for i in range(len(files)):
    xdat.append([])
    ydat.append([])
    if cdat!=None: cdat.append([])
    fname = files[i]
    with open(fname, 'r') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if row[0][0] != '#': tmp.append(row)
        xdat[i] = [x[0] for x in tmp]
        ydat[i] = [x[1] for x in tmp]
        if cdat!=None: cdat[i] = [float(x[2]) for x in tmp]

wh = (7,5)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)

gsize = None     # +gs (80,160) .. when the hexagons are not regular but taller 
# gs: (gsize)

cmap = None # eg +cmap Greys OrRd YlOrRd Oranges
# cmap: 'cmap'

clog = None
# clog: clog

cbar = None
# cbar: cbar

clabel = None
# clb: 'clabel'

for i in range(len(files)):
    args = {}
    if clog!=None: args['norm'] = mpl.colors.LogNorm()      # wrong clog labels: args['bins'] = 'log'; cb.set_label('log10(N)')
    # This version of colorbar scale requires the label
    #   if clog!=None: args['bins'] = 'log'
    #   cb.set_label('log10(N)')
    #
    # This version of colorbar does draw the dense ticks:
    if clog!=None: args['norm'] = mpl.colors.LogNorm()

    if cdat!=None: args['C'] = cdat[i]
    if cmap!=None: args['cmap'] = cmap
    if gsize!=None: args['gridsize'] = gsize
    img = ax1.hexbin(xdat[i],ydat[i],**args)
    if cbar != None:
        cb = fig.colorbar(img)
        if clabel != None: 
            cb.set_label(clabel)

xsci = None
ysci = None
# xsci: xsci
# ysci: ysci
if xsci!=None: ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')
if ysci!=None: ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')

grid = None
# gr: grid
if grid!=None: ax1.grid(color='gray', linewidth=1)

xlim = None
# xr: xlim
if xlim!=None: ax1.set_xlim(xlim)

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

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE
plt.close()

