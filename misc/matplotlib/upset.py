#!/usr/bin/env python
#
#   cat cnt.txt | mplot upset -o upset.png
#
#       1100  0.388889
#       1010  0.478261
#       0110  0.608696
#       1110  0.304348
#       1001  0.222222
#       0101  0.258065
#       1101  0.083333
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
import itertools
import csv,re,sys
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

files  = []
# FILES

fname = files[0]

labels = ''     # +lb 'lbA;lbB;lbC'
# lb: 'labels'
lb = labels.split(';')

sty = 'mine' # +sty ggplot
# sty: 'sty'

wh = (7,5)
# wh: wh

wd = 0.8
# wd: wd

# -- Read the data --
nfiles = 0
xdat = []
ydat = []
sets = {'activex':[],'activey':[],'passivex':[],'passivey':[],'lines':[]}
with open(fname) as f:
    reader = csv.reader(f, 'tab')
    tmp = []
    for row in reader:
        if row[0][0] == '#': continue
        xpos = len(xdat) - wd*0.5
        bits = list(row[0])
        nfiles = len(bits)
        min_1i = None
        max_1i = None
        for i in range(nfiles):
            if bits[i]=='1': 
                sets['activex'].append(xpos)
                sets['activey'].append(i)
                if min_1i==None: min_1i = i
                max_1i = i
            else:
                sets['passivex'].append(xpos)
                sets['passivey'].append(i)
        sets['lines'].append([[xpos,xpos],[min_1i,max_1i]])
        xdat.append(xpos)
        ydat.append(float(row[1]))

if len(lb)>0 and len(lb)!=nfiles:
    print >> sys.stderr, 'Different number of labels and sets!',lb,len(lb),nfiles
    sys.exit(1)

fcolor = None
# fc: 'fcolor'

ecolor = None
# ec: 'ecolor'

colors = [ '#337ab7', '#f0ad4e', '#5cb85c', '#5bc0de', '#d9534f', 'grey', 'black' ]

plt_args = {}       # for example: +pa "mec='grey',mfc='grey'"
# pa: {plt_args}

lbl_args = {'rotation':35,'ha':'right','multialignment':'center'} # for example: "rotation=35,ha='right',ma='center',fontsize=9"
# la: {lbl_args}


# -- Plot layout --
gs = GridSpec(2, 1, height_ratios=[5, 2])
adjust = {'left':0.12, 'right':0.95, 'bottom':0.05, 'wspace':0.05, 'hspace':0.07}
# adj: {adjust}
if adjust!=None: gs.update(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95,wspace=0.05,hspace=0.02
ax_bars = plt.subplot(gs[0])
ax_sets = plt.subplot(gs[1])

active_color  = '#444444'
passive_color = '#dddddd'

color = None
# color: 'color'
if color!=None:
    fcolor = color
    ecolor = color
    active_color = color

plt_args['color'] = active_color
if fcolor!=None: plt_args['color'] = fcolor
if ecolor!=None: plt_args['edgecolor'] = fcolor
else: plt_args['edgecolor'] = plt_args['color']
ax_bars.bar([x-wd*0.5 for x in xdat],ydat,wd,align='edge',**plt_args)


# -- Cosmetics, axes, etc --
xmin = min(xdat)-wd
xmax = max(xdat)+wd
ax_bars.set_xlim([xmin,xmax])
ax_sets.set_xlim([xmin,xmax])
ax_sets.set_ylim([-1,nfiles])
ax_bars.set_xticks([])
ax_sets.set_xticks([])
ax_sets.set_yticks([])

ylim = None
# yr: 'ylim'
if ylim!=None:                          # for example: +yr 0,1.1%
    (ymin,ymax) = ylim.split(',')
    if ymin[-1] == '%': 
        ymin = float(ymin[0:-1])*min(min(ydat))
    if ymax[-1] == '%': 
        ymax = float(ymax[0:-1])*float(max(max(ydat)))
    ax_bars.set_ylim(float(ymin),float(ymax))


# -- Sets --
ax_sets.plot(sets['activex'],sets['activey'],'.',color=active_color,mec=active_color,ms=20)
ax_sets.plot(sets['passivex'],sets['passivey'],'.',color=passive_color,mec=passive_color,ms=20)
try:
    xrange
except NameError:
    xrange = range
for i in xrange(0,nfiles,2):
    poly = Polygon([(xmin,i-0.5),(xmin,i+0.5),(xmax,i+0.5),(xmax,i-0.5)], facecolor='#eeeeee', edgecolor='#eeeeee')
    ax_sets.add_patch(poly)
for line in sets['lines']:
    ax_sets.plot(line[0],line[1],color=active_color,lw=3)
if len(lb)==nfiles:
    ax_sets.set_yticks(range(len(lb)))
    ax_sets.set_yticklabels(lb)
    ax_sets.yaxis.set_tick_params(length=0,pad=15)
    

if sty=='mine':
    for ax in [ax_bars,ax_sets]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.spines['bottom'].set_color('grey')
        ax.spines['left'].set_color('grey')
        mpl.rcParams['text.color'] = '555555'
        args = {'color':'#555555'}
        ax.patch.set_visible(False)
    for ax in [ax_sets]:
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)


xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
if xsci!=None: ax_bars.ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
if ysci!=None: ax_bars.ticklabel_format(style='sci', scilimits=ysci, axis='y')

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax_bars.set_ylabel(ylabel)
ax_bars.yaxis.labelpad = 15

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax_bars.set_xlabel(xlabel)

xticks = 1       # +xt 0
# xt: xticks
if xticks==0: ax_bars.set_xticks([])

xscale = None
# xs: 'xscale'
if xscale!=None: ax_bars.set_xscale(xscale)     # log, slog

yscale = None
# ys: 'yscale'
if yscale!=None: ax_bars.set_yscale(yscale)     # log, slog

title = None
# title: 'title'
if title!=None: ax_bars.set_title(title)

#if len(lb)>0: plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)


# dpi: dpi

# SAVE
plt.close()

