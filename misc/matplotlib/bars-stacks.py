#!/usr/bin/env python3
#
#   cat dat.txt | mplot bars-stacks -o out.png
#
#   dat:
#       wes called >500 10
#           - bar       .. if only one, each xlbl has one corresponding bar
#           - stack     .. multiple colors stacked at a single xlbl
#           - xlbl      .. placement along x-axis
#           - count     .. the y-value
#   params:
#       +o wes,acgh:called,missed_dup,missed_del                .. ordering of bars:stacks (note: xlbls order by +xt)
#       +cl -:h=//:c=#c5c5c5,ec=black:c=#f4640d:c=#007ab9       .. hatching and colors for each category defined by +o
#                                                                       .. "-" use the default
#                                                                       .. "a" alpha
#                                                                       .. "c" color, overrides ec and fc
#                                                                       .. "ec" edge color
#                                                                       .. "fc" face color
#                                                                       .. "h" hatching pattern
#                                                                       .. "hc" hatch color
#       +lg 'WES:c=...,h=...;aCGH:...'                          .. bars and stacks legend, format is the same as +cl
#       +xt '5:0-5,10:6-10,50:11-50,>500'                       .. x-ticks to display if different from xlbl data points or to enforce the order
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv,re
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

sty = None      # +sty ggplot
# sty: 'sty'
if sty=='xkcd':
    plt.xkcd()
elif sty!=None and sty!='mine':
    plt.style.use(sty)

labels = []
files  = []
# LABELS
# FILES
fname = files[0]

bars   = []
stacks = []
o_bar_stacks = None
# o: 'o_bar_stacks'
if o_bar_stacks!=None:
    tmp = o_bar_stacks.split(':')
    for x in tmp[0].split(','): bars.append(x)
    for x in tmp[1].split(','): stacks.append(x)

lb_bars   = []
lb_stacks = []
labels = None
# lb: 'labels'
if labels!=None:
    tmp = labels.split(';')
    lb_bars   = tmp[0:len(bars)]
    lb_stacks = tmp[len(bars):]


xticks = None       # +xt '0,1,2,3,4,5,10,50,100,200,500,>500'
# xt: 'xticks'
if xticks!=None:
    tmp = xticks.split(',')
    xticks = {'lbl2pos':{},'disp':[]}
    for i in range(len(tmp)):
        x = tmp[i].split(':')
        xticks['lbl2pos'][x[0]] = len(xticks['lbl2pos'])
        if len(x)>1:
            xticks['disp'].append(x[1])
        else:
            xticks['disp'].append(x[0])

dat = {}
file = open(fname,'r')
for line in file:
    row = line.rstrip('\n').split('\t')
    if row[0][0] == '#': continue
    bar   = row[0]
    stack = row[1]
    xlbl  = row[2]
    cnt   = float(row[3])
    if bar not in bars: bars.append(bar)
    if stack not in stacks: stacks.append(stack)
    if xticks==None:
        xticks={'lbl2pos':{},'disp':[]}
    if xlbl not in xticks['lbl2pos']:
        xticks['lbl2pos'][xlbl] = len(xticks['lbl2pos'])
        xticks['disp'].append(xlbl)
    xval = xticks['lbl2pos'][xlbl]      # xpos of the left bar
    if bar not in dat:
        dat[bar] = {}
    if stack not in dat[bar]:
        dat[bar][stack] = {}
        dat[bar][stack]['xval'] = []    
        dat[bar][stack]['yval'] = []
    dat[bar][stack]['xval'] +=  [0] * (xval+1-len(dat[bar][stack]['xval']))     # fill with 0's if expanding x-values to the left
    dat[bar][stack]['yval'] +=  [0] * (xval+1-len(dat[bar][stack]['yval']))
    dat[bar][stack]['xval'][xval] = xval
    dat[bar][stack]['yval'][xval] = cnt

styles = {}
for bar in dat:
    styles[bar] = {}
    for stack in dat[bar]:
        styles[bar][stack] = {}
colors = None
# cl: 'colors'
if colors!=None:
    tmp = colors.split(':')
    cl_bars   = tmp[0:len(bars)]
    cl_stacks = tmp[len(bars):]
    for i in range(len(cl_bars)):
        if cl_bars[i]=='-': continue
        bar = bars[i]
        for x in cl_bars[i].split(','):
            key,val = x.split('=')
            if key=='h':
                for stack in stacks: styles[bar][stack]['hatch'] = val
            if key=='ec':
                for stack in stacks: styles[bar][stack]['edgecolor'] = val
            if key=='hc':
                for stack in stacks: styles[bar][stack]['hatchcolor'] = val
            if key=='c':
                for stack in stacks: styles[bar][stack]['color'] = val
            if key=='a':
                for stack in stacks: styles[bar][stack]['alpha'] = float(val)
    for i in range(len(cl_stacks)):
        if cl_stacks[i]=='-': continue
        stack = stacks[i]
        for x in cl_stacks[i].split(','):
            key,val = x.split('=')
            if key=='h':
                for bar in bars: styles[bar][stack]['hatch'] = val
            if key=='hc':
                for bar in bars: styles[bar][stack]['hatchcolor'] = val
            if key=='ec':
                for bar in bars: styles[bar][stack]['edgecolor'] = val
            if key=='c':
                for bar in bars: styles[bar][stack]['color'] = val
            if key=='a':
                for bar in bars: styles[bar][stack]['alpha'] = float(val)

sty = None      # +sty ggplot
# sty: 'sty'

wh = (7,5)
# wh: wh

wd = None   # Determined automatically when not given
# wd: wd
if wd==None: wd = 1./(0.5 + len(bars))

sp = 0.1  # Small space between bars as fraction of wd; +sp 0.1
# sp: sp
if sp > 0: sp = wd*sp; wd -= sp

fig, ax1 = plt.subplots(1, 1, figsize=wh)
ytop = {}
for ibar in range(len(bars)): 
    for istack in range(len(stacks)):
        bar = bars[ibar]
        stack = stacks[istack]
        x0   = ibar*(wd+sp)
        xval = [x0+x for x in dat[bar][stack]['xval']]
        if bar not in ytop: ytop[bar] = [0.0] * len(xval)
        #print(bar,stack,dat[bar][stack],styles[bar][stack],ytop[bar])
        yval = []
        ybot = []
        for i in range(len(xval)):
            ybot.append(ytop[bar][i])
            yval.append(dat[bar][stack]['yval'][i])
            ytop[bar][i] = yval[i] + ytop[bar][i]
        args = styles[bar][stack]
        if 'edgecolor' not in args: args['linewidth'] = 0   # make optional if edgecolor set
        if 'hatchcolor' in args:
            if 'edgecolor' in args: import sys; print('Uh: not ready for both edgecolor and hatchcolor'); sys.exit(1)
            args['linewidth'] = 0
            args['edgecolor'] = args['hatchcolor']
            del(args['hatchcolor'])
        args['bottom'] = ybot
        # print(args)
        ax1.bar(xval,yval,wd,**args)


xt_args = {'rotation':35,'ha':'right','multialignment':'center'}    # +xta "rotation=35,ha='right',ma='center',fontsize=9",y=-0.05
# xta: {xt_args}

ax1.set_xticks([x+0.2 for x in range(len(xticks['disp']))])
ax1.set_xticklabels(xticks['disp'],**xt_args)

lg_args = {}            # +lga 'loc="upper left"'
# lga: {lg_args}

lgt = None
# lgt: 'lgt'
lg = None
# lg: 'lg'
if lg!=None:
    import matplotlib.patches as patches
    rct = []
    lbl = []
    for x in lg.split(';'):
        args = {}
        y = x.split(':')
        lbl.append(y[0])
        if len(y)>1:
            for x in y[1].split(','):
                z = x.split('=')
                if z[0]=='h': args['hatch'] = z[1]
                if z[0]=='ec': args['edgecolor'] = z[1]
                if z[0]=='fc': args['facecolor'] = z[1]
                if z[0]=='c': args['color'] = z[1]
                if z[0]=='a': args['alpha'] = float(z[1])
        rct.append( patches.Rectangle((1,1),1,1, **args) )
    # **lg_args does not work in < python3.5
    args = {'numpoints':1,'markerscale':2,'loc':'best','prop':{'size':10},'frameon':False,**lg_args}
    if lgt!=None: args['title'] = lgt
    plt.legend(rct,lbl,**args)

ab = None      # annotate bars:   +ab '1:aCGH;WES'
# ab: 'ab'
if ab!=None:
    tmp = ab.split(':')
    idx = int(tmp[0]) - 1
    ann = tmp[1].split(';')
    for ibar in range(len(ann)):
        xpos = idx + ibar*(wd+sp) - 0.3*wd
        ypos = ytop[bars[ibar]][idx]
        ax1.annotate(ann[ibar],xy=(xpos,ypos),ha='left',va='bottom',rotation=35)

if sty=='ggplot':
    plt.style.use('ggplot')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('grey')
    ax1.spines['left'].set_color('grey')
if sty=='mine' or sty=='xkcd':
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

hdt = None       # hide ticks: +hdt 1
# hdt: hdt
if hdt: ax1.xaxis.set_tick_params(length=0)

ylim = None
# yr: (ylim)
if ylim!=None: ax1.set_ylim(ylim)

xlim = None
# xr: (xlim)
if xlim!=None: ax1.set_xlim(xlim)

yla = {}            # ylabel args:  +yla labelpad=10
# yla: {yla}

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax1.set_ylabel(ylabel,**yla)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax1.set_xlabel(xlabel)

xscale = None
# xs: 'xscale'
if xscale!=None: ax1.set_xscale(xscale)     # log, slog

yscale = None
# ys: 'yscale'
if yscale!=None: ax1.set_yscale(yscale)     # log, slog

title = None
# title: 'title'
if title!=None: ax1.set_title(title)

#   plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)
   
adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE
plt.close()

