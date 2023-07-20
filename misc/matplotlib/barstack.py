#!/usr/bin/env python3
#
#   cat dat.txt | mplot barstack -o barplot.png +type stack-xlbl-cnt       # type/category (stacked-bar color),x-axis label,count [the default]
#
#   cat dat.txt | mplot barstack -o barplot.png +type bar-xlbl-cnt
#       dup 0   2.83
#       dup 1   30.22
#       del >1  1.78
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv,re
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

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

type = 'stack-xlbl-cnt'
# type: 'type'

direction = 'vertical'

barstacks   = []        # defines the order of the stacks
stack_label = None      # +slb 'stack1:Deletions;stack2:Duplications'           # stack labels
# slb: 'stack_label'
if stack_label!=None:
    tmp = stack_label
    stack_label = []
    barstacks   = []
    for x in tmp.split(';'):
        y = x.split(':')
        if len(y)!=2: 
            import sys
            print('Could not parse "%s", expected colon-delimited stack:label\n' % x)
            sys.exit(1)
        st = y[0]
        lb = y[1]
        barstacks.append(st)
        stack_label.append(lb)

lb = None
# xlb: 'lb'          # +xlb 'Single gene;Multiple genes'        # x-axis labels
if lb!=None: lb = lb.split(';')

xlbl2xval = {}
xticks = None       # +xt '0,1,2,3,4,5,10,50,100,200,500,>500'
# xt: 'xticks'
if xticks!=None:
    tmp = xticks.split(',')
    for i in range(len(tmp)): xlbl2xval[tmp[i]] = i

dat = {}
max_len = 0
for i in range(len(files)):
    file = open(files[i],'r')
    for line in file:
        row = line.rstrip('\n').split('\t')
        if row[0][0] == '#': continue
        stack = row[0]
        xlbl  = row[1]
        cnt   = float(row[2])
        if xlbl not in xlbl2xval:
            xlbl2xval[xlbl] = len(xlbl2xval)
        xval = xlbl2xval[xlbl]
        if stack not in dat:
            dat[stack] = {}
            dat[stack]['xval'] = []
            dat[stack]['yval'] = []
        dat[stack]['xval'] +=  [0] * (xval+1-len(dat[stack]['xval']))
        dat[stack]['yval'] +=  [0] * (xval+1-len(dat[stack]['yval']))
        dat[stack]['xval'][xval] = xval
        dat[stack]['yval'][xval] = cnt
        if stack not in barstacks: barstacks.append(stack)
        if max_len < len(dat[stack]['xval']): max_len = len(dat[stack]['xval'])

if lb==None:
    xlbl = [''] * max_len
    for lbl in xlbl2xval: xlbl[xlbl2xval[lbl]] = lbl
else:
    xlbl = lb

norm = None     # +norm frac
# norm: 'norm'
if norm=='frac' and type=='bar-xlbl-cnt':
    for stack in dat:
        sum = 0
        for y in dat[stack]['yval']: sum += y
        for j in range(len(dat[stack]['yval'])):
            dat[stack]['yval'][j] /= sum

for stack in dat:
    dat[stack]['xval'] = range(max_len)
    dat[stack]['xlbl'] = xlbl
    dat[stack]['yval'] +=  [0] * (max_len-len(dat[stack]['yval']))

sty = None      # +sty ggplot
# sty: 'sty'

wh = (7,5)
# wh: wh

wd = None   # Determined automatically when not given
# wd: wd
if wd==None: wd = 1./(2 + len(barstacks))

sp = 0  # Small space between bars as fraction of wd; +sp 0.1
# sp: sp
if sp > 0: sp = wd*sp; wd -= sp

fcolor = None
# fc: 'fcolor'

ecolor = None
# ec: 'ecolor'

# Alternative colors can be given as pre-defined colors or explicitly
#   +cl 'default'   (same as +cl '#337ab7,#f0ad4e,#5cb85c,#5bc0de,#d9534f,grey,black')
#   +cl 'sanger'
#   +cl '#01579B,#FD8230,#1B5E20,#039BE5,#9C2222,grey,black'
#
colors = None
# cl: 'colors'
if colors==None or colors=='default': colors = '#337ab7,#f0ad4e,#5cb85c,#5bc0de,#d9534f,grey,black'
elif colors=='sanger': colors = '#01579B,#FD8230,#1B5E20,#039BE5,#9C2222,grey,black'
colors = colors.split(',')

plt_args = {}       # for example: +pa "mec='grey',mfc='grey'"
# pa: {plt_args}

if direction=='vertical':
    lbl_args = {'rotation':35,'ha':'right','multialignment':'center'} # for example: "rotation=35,ha='right',ma='center',fontsize=9",y=-0.05
else:
    lbl_args = {'ha':'right','multialignment':'center'} # for example: "rotation=35,ha='right',ma='center',fontsize=9"
# la: {lbl_args}

fig, ax1 = plt.subplots(1, 1, figsize=wh)
for i in range(len(barstacks)): 
    bardat = dat[barstacks[i]]
    plt_args['color'] = colors[i%len(colors)]

    if ecolor!=None:
        plt_args['edgecolor'] = ecolor
    else:
        plt_args['linewidth'] = 0
    if fcolor!=None: plt_args['color'] = fcolor
    if stack_label!=None:
        plt_args['label'] = stack_label[i]
    else:
        plt_args['label'] = barstacks[i]

    if type=='bar-xlbl-cnt':
        xdat = [x+i*(sp+wd)-len(barstacks)*(sp+wd)*0.5 for x in bardat['xval']]
    else:
        # stack-xlbl-cnt
        xdat = [x-wd*0.5 for x in bardat['xval']]
        if i==0:
            bardat['ytop'] = bardat['yval']
        else:
            bardat['bottom'] = dat[barstacks[i-1]]['ytop']
            bardat['ytop']   = []
            for j in range(len(bardat['bottom'])):
                bardat['ytop'].append( bardat['bottom'][j] + bardat['yval'][j] )
            plt_args['bottom'] = bardat['bottom']
    ydat = bardat['yval']
    ax1.bar(xdat,ydat,wd,**plt_args)

ax1.set_xticks([x-wd*0.5 for x in dat[barstacks[0]]['xval']])
ax1.set_xticklabels(dat[barstacks[0]]['xlbl'],**lbl_args)
#ax1.xaxis.set_tick_params(length=0)

ann = None      # add annotations:   +ann '"txt":"txt","x":1,"y":155,"va":"bottom","ha":"left";"txt":"text2","x":3,"y":175'
# ann: 'ann'
if ann!=None:
    tmp = ann.split(';')
    for x in tmp:
        xdat = eval("{%s}"%x)
        txt  = xdat['txt']; del(xdat['txt'])
        xpos = xdat['x']; del(xdat['x'])
        ypos = xdat['y']; del(xdat['y'])
        ax1.annotate(txt,xy=(xpos-1,ypos),**xdat)

bh = None      # bar heights:   +ab '"va":"bottom","ha":"center","fontsize":8,"type":"int"'  .. works only with bar-xlbl-cnt atm
# bh: 'bh'
if bh!=None:
    xdat = eval("{%s}"%bh)
    intval = False
    if 'type' in xdat and xdat['type']=='int':
        intval = True
        del(xdat['type'])
    for i in range(len(barstacks)):
        bardat = dat[barstacks[i]]
        for j in range(len(bardat['yval'])):
            txt = bardat['yval'][j]
            if intval: txt = int(bardat['yval'][j])
            ax1.annotate(txt,xy=(bardat['xval'][j]+i*(sp+wd)-len(barstacks)*(sp+wd)*0.5,bardat['yval'][j]),**xdat)

if sty=='ggplot':
    plt.style.use('ggplot')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('grey')
    ax1.spines['left'].set_color('grey')
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

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax1.set_ylabel(ylabel)

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

plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)

lg_args = {}            # +lga 'loc="upper left"'
# lga: {lg_args}
lg_args = dict({'numpoints':1,'markerscale':1,'loc':'best','prop':{'size':10},'frameon':False},**lg_args)

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi

# SAVE
plt.close()

