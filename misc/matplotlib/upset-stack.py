#!/usr/bin/env python
#
#   cat venn-cnt.txt | mplot upset-stack +type venn-cnt
#       1100  0.388889
#       1010  0.478261
#
#   cat venn-stack.txt | mplot upset-stack +type venn-stack +lb 'aCGH;v3;v5'
#       1011  cds     
#       1011  cds     
#       1011  non-cds 
#
#   cat venn-bar.txt | mplot upset-stack +type venn-bar +lb 'aCGH;v3;v5' +stack 'cds:Coding;non-coding:Non-coding'     # same as venn-stack but multiple bars instead of stacked bars
#       1011  cds     
#       1011  cds     
#       1011  non-cds 
#
#   cat bar-list.txt | mplot upset-stack +type bar-list +lb 'canoes:Canoes;clamms:CLAMMS;convex:CONVEX;xhmm:XHMM' +stack 'dup:Duplications;del:Deletions'
#       dup convex,clamms     
#       dup xhmm,convex,clamms
#       del canoes

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
import itertools
import sys

type = 'venn-bar'
# type: 'type'
if type!='venn-bar' and type!='bar-list':
    print 'Only +type venn-bar and bar-list are supported atm'
    sys.exit(1)

files  = []
# FILES

fname = files[0]

labels = ''     # +lb 'lbA;lbB;lbC'
# lb: 'labels'
lb = sorted(labels.split(';'))[::-1]
list_hash = {'_list':lb}

def generate_venn_list(nlist):
    list  = []
    nvenn = (1<<nlist) - 1
    for mask in range(nvenn):
        bin_key = []
        for bit in range(nlist):
            if (1<<bit)&(mask+1):
                bin_key.append('1')
            else:
                bin_key.append('0')
        bin_key_str = ''.join(bin_key)
        list.append(bin_key_str)
    return list
    
def list2venn(list_hash,list):
    if '_nlist' not in list_hash:
        nlist = len(list_hash['_list'])
        nvenn = (1<<nlist) - 1
        for mask in range(nvenn+1):
            bin_key = []
            str_key = []
            for bit in range(nlist):
                if (1<<bit)&mask:
                    str_key.append(list_hash['_list'][bit])
                    bin_key.append('1')
                else:
                    bin_key.append('0')
            bin_key_str = ''.join(bin_key)
            str_key_str = ','.join(str_key)
            list_hash[str_key_str] = bin_key_str
    x = ','.join(sorted(list.split(','))[::-1])
    if x not in list_hash:
        nlist = len(list_hash['_list'])
        list_hash[x] = ''.join(['0'] * nlist)
        print 'Warning: the list "%s" is not known, inserting as "%s"' % (x,list_hash[x])
    return list_hash[x]


sty = 'mine' # +sty ggplot
# sty: 'sty'

wh = (5,2)  # height ratios
# wh: wh

wd = None       # determined automatically if not given
# wd: wd

# -- Read the data --
nfiles = 0
ydat = {}
udat = {}
stacks = {}
for i in range(len(files)):
    file = open(files[i],'r')
    for line in file:
        row = line.rstrip('\n').split('\t')
        if row[0][0] == '#': continue
        venn  = row[0]
        stack = 'dflt'
        value = 1.0
        if type=='bar-list':
            venn = list2venn(list_hash,row[1])
            stack = row[0]
        elif type=='venn-stack' or type=='venn-bar': stack = row[1]
        else: value = float(row[1])
        stacks[stack] = 1
        if venn in ydat:
            if stack not in ydat[venn]: ydat[venn][stack] = 0
            ydat[venn][stack] += value
            continue
        udat[venn] = {'active':[],'passive':[],'lines':[]}
        ydat[venn] = {}
        ydat[venn][stack] = value
        bits = list(venn)
        if nfiles>0 and nfiles!=len(bits):
            print >> sys.stderr, 'Inconsistent number of bits: ',line
            sys.exit(1)
        nfiles = len(bits)
        min_1i = None
        max_1i = None
        for i in range(nfiles):
            if bits[i]=='1': 
                udat[venn]['active'].append(i)
                if min_1i==None: min_1i = i
                max_1i = i
            else:
                udat[venn]['passive'].append(i)
        udat[venn]['lines'].append([min_1i,max_1i])

# plot missing:
#   - 0: leave out missing venn combinations
#   - 1: plot also missing venn combinations (implied by +cm, cumulative venn)
plot_missing_venn = 0
# pm: plot_missing_venn
if plot_missing_venn:
    venn_list = generate_venn_list(len(lb))
    for venn in venn_list:
        if venn in ydat:
            for stack in stacks:
                if stack not in ydat[venn]: ydat[venn][stack] = 0
            continue
        udat[venn] = {'active':[],'passive':[],'lines':[]}
        ydat[venn] = {}
        for stack in stacks: ydat[venn][stack] = 0
        bits = list(venn)
        nfiles = len(bits)
        min_1i = None
        max_1i = None
        for i in range(nfiles):
            if bits[i]=='1':
                udat[venn]['active'].append(i)
                if min_1i==None: min_1i = i
                max_1i = i
            else:
                udat[venn]['passive'].append(i)
        udat[venn]['lines'].append([min_1i,max_1i])

# percentage:
#   - 0: show counts
#   - 1: show percentage
#
percentage = None   # +pct 1
# pct: percentage
if percentage:
    cnt = {}
    for venn in ydat:
        for stack in ydat[venn]:
            if stack not in cnt: cnt[stack] = 0
            cnt[stack] += ydat[venn][stack]*0.01
    for venn in ydat:
        for stack in ydat[venn]:
            if cnt[stack]==0: continue
            ydat[venn][stack] /= cnt[stack]

# cumulative venn:
#   - 0:  proper Venn, disjoint sets, "10" bin does not include counts from "11"
#   - 1:  "10" bin includes counts from both "10" and "11"
#
cumulative = None   # +cm 1
# cm: cumulative
if cumulative:
    cnt = {}
    for venn in ydat:
        for venn2 in ydat:
            match = None
            for bit in range(len(venn)):
                if venn[bit]=='0': continue
                if match==None: match = 1;
                if venn2[bit]=='0': match = 0
            if match!=1: continue
            if venn not in cnt: cnt[venn] = {}
            for stack in ydat[venn]:
                if stack not in cnt[venn]: cnt[venn][stack] = 0
                cnt[venn][stack] += ydat[venn2][stack]
    for venn in cnt:
        for stack in cnt[venn]:
            ydat[venn][stack] = cnt[venn][stack]

# sort order:
#   - y: by the height of the bars [default]
#   - n: by the number of sets in the list
#
sort_order = 'y'        # +so n
# so: 'sort_order'

bars = []
for venn in ydat:
    # count_bits
    nset = 0
    for i in range(len(venn)):
        if venn[i]=='1': nset += 1

    # count y vals
    cnt  = 0
    for stack in ydat[venn]:
        cnt += ydat[venn][stack]
    if sort_order=='y':
        bars.append([venn,cnt])
    elif sort_order=='n':
        bars.append([venn,nset])
    else:
        print 'The key is not recognised: +so %s' % sort_order
        sys.exit(1)

lbl_stacks = None

stack_str = None     # +stack 'CDS:Coding;non-coding:Non-coding'
# stack: 'stack_str'

if stack_str!=None:
    stacks = []
    lbl_stacks = []
    tmp = stack_str.split(';')
    for x in tmp:
        key,val = x.split(':')
        stacks.append(key)
        lbl_stacks.append(val)
else:
    stacks = stacks.keys()

bars = sorted(bars, key=lambda x:x[1])[::-1]
for i in range(len(bars)): bars[i].append(i)

if wd==None: wd = 1./(2 + len(stacks))
dat = []
for i in range(len(stacks)):
    xval = []
    yval = []
    for bar in bars:
        xval.append(bar[2]+i*wd-len(stacks)*wd*0.5)
        venn  = bar[0]
        stack = stacks[i]
        if stack not in ydat[venn]: ydat[venn][stack] = 0
        yval.append(ydat[venn][stack])
    dat.append({'xval':xval,'yval':yval})

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
gs = GridSpec(2, 1, height_ratios=wh)
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

# plt_args['color'] = active_color
# if fcolor!=None: plt_args['color'] = fcolor
# if ecolor!=None: plt_args['edgecolor'] = fcolor
# else: plt_args['edgecolor'] = plt_args['color']
# plt_args['linewidth'] = 0

for i in range(len(dat)):
    stack = dat[i]
    plt_args['edgecolor'] = colors[i%len(colors)]
    plt_args['color']     = colors[i%len(colors)]
    plt_args['linewidth'] = 0
    if lbl_stacks!=None: plt_args['label'] = lbl_stacks[i]
    ax_bars.bar(stack['xval'],stack['yval'],wd,align='edge',**plt_args)


# -- Cosmetics, axes, etc --
xmin = 0 - len(stacks)*wd
xmax = len(bars)
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

for bar in bars:
    venn = bar[0]
    xpos = bar[2]
    ax_sets.plot([xpos for y in udat[venn]['active']],udat[venn]['active'],'.',color=active_color,mec=active_color,ms=20)
    ax_sets.plot([xpos for y in udat[venn]['passive']],udat[venn]['passive'],'.',color=passive_color,mec=passive_color,ms=20)
    for line in udat[venn]['lines']:
        ax_sets.plot([xpos for y in line],line,color=active_color,lw=3)

try:
    xrange
except NameError:
    xrange = range
for i in xrange(0,nfiles,2):
    poly = Polygon([(xmin,i-0.5),(xmin,i+0.5),(xmax,i+0.5),(xmax,i-0.5)], facecolor='#eeeeee', edgecolor='#eeeeee')
    ax_sets.add_patch(poly)
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

if lbl_stacks!=None:
    ax_bars.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)


# dpi: dpi

# SAVE
plt.close()

