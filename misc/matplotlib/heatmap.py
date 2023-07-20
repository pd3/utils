#!/usr/bin/python
#
# - not great, slow, unfinished
#
#   plot 2D distribution as heatmap
#   x y
#

# CMDLINE

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools, csv
import matplotlib.cm as cm
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

def data_coord2view_coord(p, vlen, pmin, pmax):
    dp = pmax - pmin
    dv = (p - pmin) / dp * vlen
    return dv


def nearest_neighbours(xs, ys, reso, n_neighbours):
    im = np.zeros([reso, reso])
    extent = [np.min(xs), np.max(xs), np.min(ys), np.max(ys)]
    xv = data_coord2view_coord(xs, reso, extent[0], extent[1])
    yv = data_coord2view_coord(ys, reso, extent[2], extent[3])
    for x in range(reso):
        for y in range(reso):
            xp = (xv - x)
            yp = (yv - y)
            d = np.sqrt(xp**2 + yp**2)
            im[y][x] = 1 / np.sum(d[np.argpartition(d.ravel(), n_neighbours)[:n_neighbours]])
    return im, extent

style = None        # xkcd, ggplot, ...
# sty: 'style'
if style=='xkcd':
    plt.xkcd()
elif style!=None and style!='mine':
    plt.style.use(style)

files  = []
# FILES
file = files[0]

xdat = []
ydat = []
with open(file, 'rb') as f:
    reader = csv.reader(f, 'tab')
    for row in reader:
        if row[0][0] != '#':
            xdat.append(float(row[0]))
            ydat.append(float(row[1]))

wh = (7,5)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)

cmap = None # eg +cmap Greys
# cmap: 'cmap'

clog = 0
# clog: clog

cbar = 0
# cbar: cbar

resolution = 256
# rs: resolution

neighbours = 32
# nb: neighbours

args = {}
if clog>0: args['bins'] = 'log'
if cmap!=None: args['cmap'] = cmap
im, extent = nearest_neighbours(xdat, ydat, resolution, neighbours)
ax1.imshow(im, origin='lower') #, extent=extent, cmap=cm.jet)


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

