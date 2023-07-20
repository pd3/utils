#!/usr/bin/env python3
#
#   cat chr-x-y.txt         | mplot manhattan -o img.png
#   cat chr-x-y-color.txt   | mplot manhattan -o img.png +type chr-x-y-col      # "." for default
#
# CMDLINE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math

style = None        # xkcd, ggplot, ...
# sty: 'style'
if style=='xkcd':
    plt.xkcd()
elif style!=None and style!='mine':
    plt.style.use(style)

files  = []
# FILES
file = files[0]

type = 'chr-x-y'
# type: 'type'

th_args = dict(color='#f4640d')
th = None       # +th '0.05,linestyle="--",lw=1'
# th: 'th'
if th!=None:
    tmp = th.split(',')
    th = float(tmp[0])
    if len(tmp)>1:
        xx = eval("dict(%s)"%','.join(tmp[1:]))
        th_args = { **th_args, **xx }

# Alternative colors can be given as pre-defined colors or explicitly
#   +cl 'default'   (same as +cl '#337ab7,#f0ad4e,#5cb85c,#5bc0de,#d9534f,grey,black')
#   +cl 'sanger'
#   +cl '#01579B,#FD8230,#1B5E20,#039BE5,#9C2222,grey,black'
#   +cl grey
#
colors = None
# cl: 'colors'
if colors==None or colors=='default': colors = '#007ab9,#faa300,#00b7ed,#242424,#00a077,#f5e636,#f4640d,#e37dac'       # colors = '#337ab7,#f0ad4e,#5cb85c,#5bc0de,#d9534f,grey,black'
elif colors=='grey': colors = '#474747,#b5b5b5'
elif colors=='sanger': colors = '#01579B,#FD8230,#1B5E20,#039BE5,#9C2222,grey,black'
colors = colors.split(',')

try:
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)
except:
    # deprecated: 
    mpl.rcParams['axes.color_cycle'] = colors

def chrpos2x(dat,chr,pos,color=None):
    global colors
    if 'off' not in dat:
        dat['off']  = {chr:0}
        dat['prev_chr'] = chr
        dat['prev_pos'] = pos
        dat['chrpos']   = 0
        dat['colors']   = []
        dat['xticks']   = [0]
        dat['xlbls']    = [chr]
    if dat['prev_chr']!=chr:
        if chr in dat['off']: 
            print('Chromosomes not in blocks, see e.g. '+chr)
            sys.exit(1)
        dat['off'][chr] = dat['chrpos']
        dat['prev_pos'] = 0
        dat['xticks'].append(dat['chrpos'])
        dat['xlbls'].append(chr)
    if pos < dat['prev_pos']:
        print('The file is not sorted, see e.g. '+chr+':'+pos)
        sys.exit(1)
    if color!=None:
        dat['colors'].append(color)
    elif len(dat['off'])%2:
        dat['colors'].append(colors[0])
    else:
        dat['colors'].append(colors[1])
    dat['xticks'][-1] = dat['off'][chr] + 0.5*pos
    dat['prev_pos'] = pos
    dat['prev_chr'] = chr
    dat['chrpos'] = dat['off'][chr] + pos
    return dat['chrpos']

xdat = []
ydat = []
chrs = {}
f = open(file, 'r')
for line in f:
    row = line.split("\t")
    row[-1] = row[-1].rstrip('\n')
    pval = float(row[2])
    color = None
    if type=='chr-x-y' and th!=None and pval<=th: color = th_args['color']      # color points above the line in the chr-x-y mode
    if type=='chr-x-y-col' and row[3]!='.': color = row[3]                      # color provided explicitly by the caller
    xval = chrpos2x(chrs,row[0],float(row[1]),color)
    xdat.append(xval)
    ydat.append(pval)  # pvalue

# sort the dots so that the colored are on top
sort_colors = 1
if sort_colors==1:
    from functools import cmp_to_key
    cdat = chrs['colors']
    idx = []
    for i in range(len(xdat)):
        idx.append(i)
    def cmp_colors(a, b):
        global cdat,colors
        if cdat[a]==cdat[b]: return 0
        adflt = 0
        bdflt = 0
        if cdat[a]==colors[0] or cdat[a]==colors[1]: adflt = 1
        if cdat[b]==colors[0] or cdat[b]==colors[1]: bdflt = 1
        if adflt==1 and bdflt==1: return 0
        if adflt < bdflt: return 1
        if adflt > bdflt: return -1
        return 0
    sidx = sorted(idx, key=cmp_to_key(cmp_colors))
    xdat_new = []
    ydat_new = []
    cdat_new = []
    for i in range(len(sidx)):
        xdat_new.append(xdat[sidx[i]])
        ydat_new.append(ydat[sidx[i]])
        cdat_new.append(cdat[sidx[i]])
    xdat = xdat_new
    ydat = ydat_new
    chrs['colors'] = cdat_new

yscale = None       # +ys -log10,log,symlog
# ys: 'yscale'

wh = (7,3)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)

plt_args = {}       # for example: +pa "mec='grey',mfc='grey',zorder=100,clip_on=False,alpha=0.5"
if style=='mine': plt_args = {'zorder':100,'clip_on':False,'s':10}
# pa: {plt_args}

if yscale!=None:
    if yscale=='-log10':
        for i in range(len(ydat)): ydat[i] = -math.log10(ydat[i])
    else:
        ax1.set_yscale(yscale)

ax1.scatter(xdat,ydat,c=chrs['colors'],**plt_args)

if th!=None:
    if yscale=='-log10':
        lth = -math.log10(th)
    th_args['color'] = colors[1]    # grey dashed line
    ax1.plot([min(xdat),max(xdat)],[lth,lth],**th_args)

ysci = None         # +ysci -2,2
# ysci: (ysci)
if ysci!=None: ax1.ticklabel_format(style='sci', scilimits=ysci, axis='y')

ax1.set_xticks(chrs['xticks'])
ax1.set_xticklabels(chrs['xlbls'])
ax1.xaxis.set_tick_params(length=0)

# remove overlapping labels
plt.gcf().canvas.draw()
ticks = [tick for tick in plt.gca().get_xticklabels()]
bnd = 0
for i, t in enumerate(ticks):
    bb = t.get_window_extent()
    if bb.x0 < bnd:
        t.update(dict(alpha=0))
    else:
        bnd = bb.x1 + (bb.x1-bb.x0)*0.5
    #print("Label ", i, ", data: ", t.get_text(), " ; ", t.get_window_extent())

ylabel = None       # +yl label, +yl label:4
# yl: 'ylabel'
if ylabel!=None:
    pad = 4
    import re
    m = re.search('(.+):(\d+)',ylabel)
    if m!=None:
        ylabel = m.group(1)
        pad = float(m.group(2))
    ax1.set_ylabel(ylabel,labelpad=pad)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel,labelpad=10)

ylim = None 
# yr: 'ylim'
if ylim!=None:                          # for example: +yr 0,1.1%
    (ymin,ymax) = ylim.split(',')
    if ymin[-1] == '%':
        ymin = float(ymin[0:-1])*min(min(ydat))
    if ymax[-1] == '%':
        ymax = float(ymax[0:-1])*float(max(max(ydat)))
    ax1.set_ylim(float(ymin),float(ymax))


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

#   lg_args = {}
#   # lga: {lg_args}
#   lg_args = dict({'numpoints':1,'markerscale':1,'loc':'best','prop':{'size':10},'frameon':False},**lg_args)
#   
#   if len(lb)>0: plt.legend(**lg_args)
#   

adjust = dict(left=0.08,right=0.99,bottom=0.15)
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95


# dpi: dpi

# SAVE
plt.close()

