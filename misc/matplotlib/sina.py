#!/usr/bin/env python3
#
#   # +type xyc
#   #   - x are strings correspondig to centers of the violins (as if individual bars in a barplot)
#   #   - y is the y-axis placement
#   #   - c is [0-1] value indicating the shade
#   cat xyc.txt | mplot sina +type xyc
#       snv   10    0.5
#       snv   100   1
#       snv   20    0
#       indel 20    0.3
#       indel 5     0.1
#
#   cat xy.txt | mplot sina +type xy        # same as +type xyc but using 1 for shading
#   cat xy.txt | mplot sina +type xycb      # same as +type xyc but overlay a boxplot created from 4-th column
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
import matplotlib.cm as cm
import sys,random

fs = None   # +fs 18    # increase font size of all texts
# fs: fs
if fs!=None: plt.rcParams.update({'font.size': fs})

type = 'xyc'    # +type <type>  .. xyc xy xycb
# type: 'type'
if type!='xyc' and type!='xy' and type!='xycb':
    print('Only +type xyc, xy, or xycb is supported atm')
    sys.exit(1)

do_color = False
do_box   = False
if type=='xycb': do_box = True
if type=='xycb' or type=='xyc': do_color = True

files  = []
# FILES
fname = files[0]

sty = 'mine' # +sty ggplot
# sty: 'sty'

wh = (5,4)  # height ratios
# wh: wh

nbin = 100     # +n 100
# n: nbin

jitter = ['0,0']  # +jr 0.5,0.5       # x=fraction of the maximum value,y=abs,y%=rel
# jr: 'jitter'
jitter = jitter.split(',')

xmin = None
xmax = None

def color_scale(c1,c2,mix=0):
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def percentile(vals,p):
    N = len(vals)
    n = p*(N+1)
    k = int(n)
    d = n-k
    if k<=0: return vals[0]
    if k>=N: return vals[N-1]
    return vals[k-1] + d*(vals[k] - vals[k-1])

def adjacent_values(vals):
    q1 = percentile(sdat,0.25)
    q3 = percentile(sdat,0.75)
    uav = q3 + (q3-q1)*1.5
    if uav > vals[-1]: uav = vals[-1]
    if uav < q3: uav = q3
    lav = q1 - (q3-q1)*1.5
    if lav < vals[0]: lav = vals[0]
    if lav > q1: lav = q1
    return [lav,uav]

# -- Read the data --
def read_data(fname):
    global xmin,xmax,jitter
    ymin = ymax = None
    dat  = { 'raw':{}, 'binned':{}, 'box':{}, 'xticks':[], 'xtick_labels':[] }
    file = open(fname,'r')
    for line in file:
        row = line.rstrip('\n').split('\t')
        if row[0][0] == '#': continue
        bar = row[0]
        if bar not in dat['raw']:
            dat['raw'][bar] = []
            dat['box'][bar] = []
            dat['xticks'].append(len(dat['raw'])+1)
            dat['xtick_labels'].append(row[0])
        y = float(row[1])
        if ymin==None or ymin>y: ymin = y
        if ymax==None or ymax<y: ymax = y
        c = 1
        if do_color: c = float(row[2])
        dat['raw'][bar].append({'y':y,'c':c})
        if do_box: dat['box'][bar].append(float(row[3]))

    for bar in dat['raw']:
        bins = [{'y':0,'n':0,'c':[]} for x in range(nbin+1)]
        raw  = dat['raw'][bar]
        for row in raw:
            i = int(nbin*(row['y']-ymin)/(ymax-ymin))
            bins[i]['y'] = row['y']
            bins[i]['n'] += 1
            bins[i]['c'].append(row['c'])
        dat['binned'][bar] = bins

    nmax = 0
    for bar in dat['binned']:
        bins = dat['binned'][bar]
        for bin in bins:
            if nmax < bin['n']: nmax = bin['n']

    xmin = 0.5*nmax
    xmax = len(dat['binned'])*nmax+0.5*nmax
    xjitter = nmax*float(jitter[0])
    yjitter = jitter[1]
    if yjitter[-1]=='%': yjitter = -float(yjitter[:-1])
    else: yjitter = float(yjitter)

    for i in range(len(dat['xticks'])):
        dat['xticks'][i] = (dat['xticks'][i]-1)*nmax

    xdat = []
    ydat = []
    cdat = []
    bdat = {'dat':[],'pos':[],'wd':[]}
    ibar = 0
    for bar in dat['binned']:
        ibar += 1
        bins = dat['binned'][bar]
        for bin in bins:
            icolor = 0
            x0 = -0.5*bin['n'] + ibar*nmax
            cvals = sorted(bin['c'],reverse=True)
            for i in range(bin['n']):
                xrand = random.random()*xjitter - 0.5*xjitter
                yjmax = yjitter
                if yjmax<0: bin['y']*(-yjmax)/100.
                yrand = random.random()*yjmax - 0.5*yjmax
                # this can only be used when colorbar is not needed, otherwise cmap has to be created:
                #   color = color_scale(colors[0],colors[1],cvals[i])
                color = cvals[i]
                xdat.append(xrand + x0+i)
                ydat.append(yrand + bin['y'])
                cdat.append(color)
        if do_box:
            bdat['dat'].append(dat['box'][bar])
            bdat['pos'].append(ibar*nmax)
            bdat['wd'].append(xmax*0.15)
    return (dat,xdat,ydat,cdat,bdat)

dat,xdat,ydat,cdat,bdat = read_data(fname)

# -- create colormap --

colors = None   # +cl '#242424,#faa300,#f4640d'     # 3rd color for boxplot, optional
# cl: 'colors'
if colors!=None:
    colors = colors.split(',')
    if len(colors)==1: colors[2] = colors[1] = colors[0]
else:
    colors = ['#242424','#faa300','#f4640d']

clo = mpl.colors.to_rgb(colors[0])
chi = mpl.colors.to_rgb(colors[1])
cm.register_cmap(cmap=mpl.colors.LinearSegmentedColormap('mycmap',{
    'red':   [(0.,clo[0],clo[0]), (1.,chi[0],chi[0])],
    'green': [(0.,clo[1],clo[1]), (1.,chi[1],chi[1])],
    'blue':  [(0.,clo[2],clo[2]), (1.,chi[2],chi[2])],
    }))

# -- plot the data --

plt_args = {}       # for example: +pa "mec='grey',mfc='grey',zorder=100,clip_on=False,alpha=0.5"
plt_args = {'zorder':100,'clip_on':False,'alpha':0.5}
# pa: {plt_args}


fig, ax1 = plt.subplots(1, 1, figsize=wh)
sc = ax1.scatter(xdat,ydat, s=10,marker='o',c=cdat,cmap='mycmap',**plt_args)
if do_box:
    parts = ax1.boxplot(bdat['dat'],positions=bdat['pos'],widths=bdat['wd'],patch_artist=True,vert=True,showfliers=False)
    ec = 'black'
    for i in range(len(parts['boxes'])):
        pc = parts['boxes'][i]
        fc = [1,1,1,0]
        ec = colors[2]
        pc.set_facecolor(fc)
        pc.set_edgecolor(ec)
    for item in ['whiskers', 'fliers', 'medians', 'caps']:
        for pc in parts[item]:
            pc.set_color(ec)
            pc.set_linewidth(1.5)


cbl = None       # +cbl 'Colorbar label'
# cbl: 'cbl'
if cbl!=None:
    # https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.colorbar.html
    cb = fig.colorbar(sc,shrink=0.7,aspect=15,pad=0.0) #, orientation='horizontal',#aspect=10, pad=0.04)
    cb.set_label(cbl,labelpad=10)
    cb.set_clim(0,1)

if sty=='mine':
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('grey')
    ax1.spines['left'].set_color('grey')
    mpl.rcParams['text.color'] = '555555'
    #args = {'color':'#555555'}
    ax1.patch.set_visible(False)

yscale = None       # +ys log,symlog
# ys: 'yscale'
if yscale!=None: ax1.set_yscale(yscale)     # log, symlog

ta_args = {'y':1.08}       # for example: +ta "y=1.08"
# ta: {ta_args}
title = None
# title: 'title'
if title!=None: ax1.set_title(title,**ta_args)

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax1.set_ylabel(ylabel,labelpad=5)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel,labelpad=10)

xlim = None     # +xr 1.5   # extend nmax (+10% in this example)
# xr: xlim
if xlim==None: xlim = 1.1
x = 0.5*(xmax+xmin)
w = (xmax-xmin)*xlim
ax1.set_xlim(x-w*0.5,x+w*0.5)

xsci = None
ysci = None
# xsci: (xsci)
# ysci: (ysci)
if xsci!=None: ax1.ticklabel_format(style='sci', scilimits=xsci, axis='x')        # +xsci -2,2
if ysci!=None: ax1.ticklabel_format(style='sci', scilimits=ysci, axis='y')

lbl_args = {'rotation':35,'ha':'right','multialignment':'center'}   # for example: +la "rotation=35,ha='right',ma='center',fontsize=9"
# la: {lbl_args}

ax1.set_xticks(dat['xticks'])
ax1.set_xticklabels(dat['xtick_labels'],**lbl_args)

hdt = None       # hide ticks: +hdt 1
# hdt: hdt
if hdt:
    ax1.xaxis.set_tick_params(length=0)


adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: +adj bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE

plt.close()

