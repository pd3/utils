#!/usr/bin/env python3
#
#   cat type-cnt.txt    | mplot tstv -o barplot.png +type xcnt
#   cat type.txt        | mplot tstv -o barplot.png +type x
#
#       T>A	1038
#       A>G	4622
#       ...
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys

style = 'mine'
# sty: 'style'

type = 'xcnt'
# type: 'type'

files  = []
# FILES

keys = {
   'A>C':0,
   'A>G':1,
   'A>T':2,
   'C>A':3,
   'C>G':4,
   'C>T':5,
   'G>A':6,
   'G>C':7,
   'G>T':8,
   'T>A':9,
   'T>C':10,
   'T>G':11
}

ts = [ 'A>G','G>A','C>T','T>C' ]
nts = 0
ntv = 0

dat = []
for key in keys: dat.append([keys[key],key,0])

with open(files[0], 'r') as f:
    tmp = []
    for line in f:
        row = line.split("\t")
        row[-1] = row[-1].rstrip('\n')
        if row[0] not in keys: sys.exit("Unknown key: ",row[0])
        idx = keys[row[0]]
        cnt = 0
        if type=='xcnt': cnt = float(row[1])
        else: cnt = 1
        dat[idx][2] += cnt
        if row[0] in ts: nts += cnt
        else: ntv += cnt

n = 12
col  = list(range(n))
ecol = list(range(n))
for i in range(n):
    col[i]  = '#d9534f'  # red
    ecol[i] = '#D43F3A'
col[1]  = col[5]  = col[6]  = col[10]  = '#337ab7' # blue
ecol[1] = ecol[5] = ecol[6] = ecol[10] = '#2E6DA4' # blue

wh = (5,4)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)

ax1.bar([row[0] for row in dat], [row[2] for row in dat], color=col, edgecolor=ecol)

ylabel = 'Count'
# yl: 'ylabel'
ax1.set_ylabel(ylabel)

title = None
# title: 'title'
if title!=None: ax1.set_title(title)
elif ntv!=0: ax1.set_title("ts/tv=%.2f, n=%d" % (nts/ntv,nts+ntv))

ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
ax1.set_xlim(-0.5,n+0.5)
plt.xticks([row[0] for row in dat],[row[1] for row in dat],rotation=45,fontsize=9)
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

