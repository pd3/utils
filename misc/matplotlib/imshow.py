#!/usr/bin/python
#
#   - X1  X2 ..          Xn           # first line is x-axis range if `+xa 1` is given
#   Y1  y11 y12 y13 y14 .. y1n
#   ..  ..
#   Ym  ym1 ym2 ym3 ym4 .. ymn
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import itertools
import csv
import math
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

dat    = []
files  = []
labels = []
# FILES
# LABELS

xa = None       # first line is x-axis range if `+xa 1` is given
# xa: xa

xlim = None
# xr: xlim

hl = []
# hl: [hl]

xs = None
# xs: xs

cmap = 'Blues' # +cmap terrain_r, Spectral_r, rdylgn_r, RdYlBu_r, Greens, gnuplot_r, coolwarm, BuGn
# cmap: 'cmap'

dat = []
ymin = []
ymax = []
iext = None
for i in range(len(files)):
    dat.append([])
    fname = files[i]
    with open(fname, 'r') as f:
        reader = csv.reader(f, 'tab')
        tmp = []
        for row in reader:
            if xa!=None and iext==None:
                iext = [1,len(row)]
                if xlim!=None:
                    for j in range(1,len(row)):
                        if float(row[j])<float(xlim[0]): iext[0] = j;
                        if float(row[j])<float(xlim[1]): iext[1] = j
                else:
                    xlim = [float(row[iext[0]]),float(row[iext[1]-1])]
                continue
            if row[0][0] == '#': continue
            if len(ymin)<=i: 
                ymin.append(float(row[0]))
                ymax.append(float(row[0]))
            if ymin[i] > row[0]: ymin[i] = float(row[0])
            if ymax[i] < row[0]: ymax[i] = float(row[0])
            tmpx = []
            max  = 0
            for j in range(iext[0],iext[1]):
                if float(row[j]) > max: max = float(row[j])
                tmpx.append(float(row[j]))
            if xs!=None and max!=0:
                for j in range(len(tmpx)): tmpx[j] = tmpx[j] / max
            tmp.append(tmpx)
        dat[i] = tmp

        # m = max(tmp)
        # if len(dmax)<=i: dmax.append(m)
        # elif m<dmax[i]: dmax[i] = m
        # m = min(tmp)
        # if len(dmin)<=i: dmin.append(m)
        # elif m<dmin[i]: dmin[i] = m

width  = 2
height = 2.3
ncols  = 2
nrows  = 1
if len(files)==1: ncols = 1; width = 7; height = 7;
# ncols
# nrows
# wd: width
# ht: height

ylabel = None
# yl: 'ylabel'

xlabel = None
# xl: 'xlabel'

labelsize = None
# ls: labelsize

colorbar = None
# cbar: colorbar

interpolation = 'none'  # +intpl nearest, bilinear, etc; default is none
# intpl: 'interpolation'

fig, ax = plt.subplots(nrows, ncols, figsize=(ncols*width,nrows*height))
for i in range(ncols*nrows):
    irow = i / ncols; icol = i % ncols
    if i>=len(labels): ax[irow,icol].axis('off'); continue;
    if len(files)==1: tax = ax
    else: tax = ax[irow,icol]
    args = {}
    args['origin'] = 'lower'
    args['interpolation'] = interpolation
    args['aspect'] = 'auto'
    if cmap!=None: args['cmap'] = cmap
    extent = [0,len(dat[i][0]),ymin[i],ymax[i]]
    if xlim!=None:
        extent[0] = xlim[0]
        extent[1] = xlim[1]
    args['extent'] = extent
    #args['norm'] = LogNorm(vmin=dmin[i],vmax=dmax[i])
    im = tax.imshow(dat[i], **args)
    if len(hl)>0: 
        fig.canvas.draw()
        box = ax._position.bounds
        tmpax = fig.add_axes([box[0], box[1], box[2], box[3]])
        tmpax.set_axis_off()
        tmpax.set_ylim(0,1)
        for line in hl:
            y = (line-ymin[i])/(ymax[i]-ymin[i])
            tmpax.plot([0,1],[y,y],'--',lw=1.8,color='red')
    if ylabel!=None: tax.set_ylabel(ylabel,fontsize=labelsize)
    if xlabel!=None: tax.set_xlabel(xlabel,fontsize=labelsize)
    if colorbar!=None: plt.colorbar(im,ax=tax)

title = None
# title: 'title'
if title!=None: ax.set_title(title)

xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
if xsci!=None: ax.ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
if ysci!=None: ax.ticklabel_format(style='sci', scilimits=ysci, axis='y')

#plt.subplots_adjust(bottom=0.05,left=0.05,right=0.95,top=0.95,hspace=0.12,wspace=0.12)

# SAVE
plt.close()

