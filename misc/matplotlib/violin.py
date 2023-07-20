#!/usr/bin/env python
#
# zcat lengths.txt.gz | mplot violin -o plots/lengths.png
# zcat lengths.txt.gz | mplot violin -o plots/lengths.png +lb 'XY:Label 1;ZY:Label 2;;WZ:Label 3, padded'
#
# - arbitrary nuber of categories, create violin distribution for each
#
# For example:
#   LINE    266
#   LINE    483
#   SINE    102
#   ..
#
# Options:
#       +ps 1       pre-sorted, use violin in the order givne
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv,sys,random
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

style = None        # xkcd, ggplot, ...
# sty: 'style'
if style=='xkcd':
    plt.xkcd()
elif style!=None and style!='mine':
    plt.style.use(style)

files = []
# FILES

order = []
groups = {}
for i in range(len(files)):
    fname = files[i]
    with open(fname) as f:
        reader = csv.reader(f, 'tab')
        for row in reader:
            if row[0][0] == '#': continue
            row[0] = row[0].replace('\\n','\n')
            if row[0] not in groups:
                groups[row[0]] = []
                order.append(row[0])
            groups[row[0]].append(float(row[1]))

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

def trim_at_percentile(vals,p):
    if p[0]>1 or p[1]>1:
        p[0] = p[0]/100.
        p[1] = p[1]/100.
    idx1 = int(len(vals)*float(p[0]))
    idx2 = int(len(vals)*float(p[1]))
    return vals[idx1:idx2]

trim_pctl = None        # +trim (1,99) or +trim (0.01,0.99)
# trim: [trim_pctl]

presorted = 0          # +ps 1
# ps: [presorted]

keys = []
if presorted==0: keys = sorted(groups.keys())
else: keys = order
lab  = keys    # labels
pos  = []

lb = None
# lb: 'lb'
if lb!=None:
    keys = []
    lab  = []
    key_lbs = lb.split(';')
    ipos = 0
    for key_lb in key_lbs:
        if key_lb=='': ipos += 1; continue
        ipos += 1
        key,lb = key_lb.split(':')
        lab.append(lb)
        keys.append(key)
        pos.append(ipos)
else:
    for key in keys: pos.append(len(pos))
    

med = []    # medians
iqr = []    # inter-quantile ranges
avs = []    # upper and lower adjacent values
dat = []
for grp in keys:
    #print grp
    sdat = sorted(groups[grp])
    if trim_pctl!=None: sdat = trim_at_percentile(sdat,trim_pctl)
    med.append(percentile(sdat,0.5))
    iqr.append([percentile(sdat,0.25),percentile(sdat,0.75)])
    avs.append(adjacent_values(sdat))
    dat.append(sdat)

wh = (7,5)
# wh: wh

plot_med = None     # +med "c='white',ms=6" or "c=\'white\',ms=6" if run through wrapper
# med: {plot_med}

plot_iqr = None     # +iqr "c='black',lw=5"     .. inter-quantile range (25-75th percentile, thick bar)
# iqr: {plot_iqr}

plot_avs = None     # +avs "c='black',lw=1"     .. adjacent values (1.5*(75-25th percentile), thin bars)
# avs: {plot_avs}

ysci = None
# ysci: (ysci)  # +ysci -2,2

fcolor = '#D43F3A'      # red; +fc 'convex:#D43F3A,clamms:#D43F3A,..'
# fc: 'fcolor'
fcolors = fcolor.split(',')
if len(fcolors)==1:
    fcolors = [fcolor] * len(keys)
else:
    tmp = {}
    for x in fcolors:
        key,col = x.split(':')
        tmp[key] = col
    fcolors = []
    for grp in keys:
        fcolors.append(tmp[grp])

ecolor = 'black'
# ec: 'ecolor'

plot_dat = 0        # show data points? +pd 1
# pd: plot_dat

fig, ax = plt.subplots(1, 1, figsize=wh)

try:
    parts = ax.violinplot(dat,pos,showmeans=False,showmedians=False,showextrema=False)
    if plot_dat:
        for i in range(len(dat)): ax.plot([i + 0.25*random.random() - 0.125 for x in dat[i]],dat[i],'.',color=ecolor,alpha=0.5)
except Exception as e:
    print(e)
    plot_med = {'c':'black','ms':6}
    parts = {}
    parts['bodies'] = []
for i in range(len(med)):
    if plot_avs!=None: ax.plot([i,i],avs[i],'-',**plot_avs)
    if plot_iqr!=None: ax.plot([i,i],iqr[i],'-',**plot_iqr)
    if plot_med!=None: ax.plot(i,med[i],'o',mec='none',**plot_med)
for i in range(len(parts['bodies'])):
    pc = parts['bodies'][i]
    fc = fcolors[i]
    pc.set_facecolor(fc)
    pc.set_edgecolor(ecolor)
    pc.set_alpha(1)

xlab_args = {}      # for example: +xarg "rotation=45,ha='right',ma='left',fontsize=9"
# xarg: {xlab_args}         
ax.get_xaxis().set_tick_params(direction='out')
ax.xaxis.set_ticks_position('bottom')
ax.set_xticks(pos)
ax.set_xticklabels(lab,**xlab_args)
ax.xaxis.set_ticks_position('none')

if ysci!=None: ax.ticklabel_format(style='sci', scilimits=ysci, axis='y')
#if len(labels): ax.set_ylabel(labels[i])

yscale = None
# ys: 'yscale'
if yscale!=None: ax.set_yscale(yscale)     # log, slog

ylim = None
# yr: ylim
if ylim!=None: ax.set_ylim(ylim)

xlabel = None
# xl: 'xlabel'
if xlabel!=None: ax.set_xlabel(xlabel)

ylabel = None
# yl: 'ylabel'
if ylabel!=None: ax.set_ylabel(ylabel)

ax.set_xlim(pos[0]-0.75,pos[-1]+0.75)

args = {}
if style=='ggplot':
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines['bottom'].set_color('grey')
    ax.spines['left'].set_color('grey')
    mpl.rcParams['text.color'] = '555555'
    args = {'color':'#555555'}
if style=='mine':
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines['bottom'].set_color('grey')
    ax.spines['left'].set_color('grey')
    mpl.rcParams['text.color'] = '555555'
    args = {'color':'#555555'}
    ax.patch.set_visible(False)

title = None
# title: 'title'
if title!=None: ax.set_title(title)

adjust = {'left':0.08,'bottom':0.08,'right':0.95,'top':0.95}
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# SAVE
plt.close()

