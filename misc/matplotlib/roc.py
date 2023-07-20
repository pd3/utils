#!/usr/bin/env python3
#
#   cat val-is_tp.txt | mplot roc -o roc.png +fn 105
#
#       1323 0
#       1324 1
#       ...
#
#   options:
#       +fn 1500
#           number of false negatives not present in the data set
#
#       +pa 'marker="o",ls="-"'
#           plot arguments, e.g.
#
#       +nth 10
#           number of threshold values to annotate the curve with
#
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys

style = 'mine'
# sty: 'style'

files  = []
# FILES

nFP_tot = 0
nTP_tot = 0
nFN_tot = 0      # +fn 105
# fn: nFN_tot

dat = []
with open(files[0], 'r') as f:
    tmp = []
    for line in f:
        row = line.split("\t")
        row[-1] = row[-1].rstrip('\n')
        value = float(row[0])
        is_tp = int(row[1])
        if is_tp!=0 and is_tp!=1: sys.exit('Expected [01] in the is_tp column, found: ',line)
        if is_tp==1: nTP_tot += 1
        else: nFP_tot += 1
        dat.append([value,is_tp])

sdat = sorted(dat, key=lambda x:x[0], reverse=True)
dat  = []
prev = None
nFP  = 0
nTP  = 0
for x in sdat:
    if x[1]==1: nTP += 1
    else: nFP += 1
    if prev==x[0]: continue
    prev = x[0]
    xval = nFP / (nFP_tot+nTP_tot)      # FP rate
    yval = nTP / (nFN_tot+nTP_tot)      # TP rate
    dat.append([xval,yval])

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

plt_args = [{}]       # for example: +pa "mec='grey',mfc='grey',zorder=100,clip_on=False,alpha=0.5"
# pa: [{plt_args}]

wh = (5,4)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)

for i in range(len(files)):
    args = plt_args[i]
    if style=='mine': args = dict({'zorder':100,'clip_on':False},**args)
    ax1.plot([row[0] for row in dat], [row[1] for row in dat],**args)

xlabel = 'FP rate'
# xl: 'xlabel'
ax1.set_xlabel(xlabel)

ylabel = 'TP rate'
# yl: 'ylabel'
ax1.set_ylabel(ylabel)

ylim = None
# yr: 'ylim'
if ylim!=None:
    x = ylim.split(',')
    if x[0]!='': ax1.set_ylim(bottom=float(x[0]))
    if x[1]!='': ax1.set_ylim(top=float(x[1]))

xlim = None         # for example: +xr 0,1.1%
# xr: 'xlim'
if xlim!=None:
    x = xlim.split(',')
    if x[1][-1] == '%':
        x[1] = float(x[1][0:-1])*max([row[0] for row in dat])
    if x[0]!='': ax1.set_xlim(left=float(x[0]))
    if x[1]!='': ax1.set_xlim(right=float(x[1]))

ta_args = {}        #  +ta t=1.0
# ta: {ta_args}
title = None
# title: 'title'
if title!=None: ax1.set_title(title,**ta_args)

plt.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.9)

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

# dpi: dpi

# SAVE
plt.close()

