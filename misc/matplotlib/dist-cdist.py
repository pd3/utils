#!/usr/bin/python
#
#   Input:
#       beg  end  cnt 
#
#   Notes:
#       [beg,end) intervals, plots density and cumulative distribution
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import csv
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

labels = []
files  = []
# LABELS
# FILES


xdat = []
ydat = []
cdat = []
for i in range(len(files)):
    xdat.append([])
    ydat.append([])
    cdat.append([])
    fname = files[i]
    with open(fname, 'rb') as f:
        reader = csv.reader(f, 'tab')
        sum = 0
        for row in reader:
            if row[0][0] == '#': continue
            beg = float(row[0])
            end = float(row[1])
            cnt = float(row[2])
            xdat[-1].append(beg)
            ydat[-1].append(cnt/(end-beg))
            sum += cnt
            cdat[-1].append(sum)
        for i in range(len(cdat[-1])):
            cdat[-1][i] = cdat[-1][i]/sum

line1  = '.-'
line2  = '.-'
label1 = 'Cumulative fraction'
label2 = 'Density'
col1 = 'black'
col2 = '#D43F3A'

fig, ax1 = plt.subplots(1, 1, figsize=(7,5))
ax2 = ax1.twinx()
plots = []
for i in range(len(labels)): 
    plots += ax1.plot(xdat[i],cdat[i],'.-',color=col1,label=label1)
    plots += ax2.plot(xdat[i],ydat[i],'.-',color=col2,label=label2)

ax2.ticklabel_format(style='sci', scilimits=(-2,2), axis='y')
ax2.ticklabel_format(style='sci', scilimits=(-2,2), axis='x')

ax1.set_ylabel(label1, color=col1)
ax2.set_ylabel(label2, color=col2)

for tl in ax1.get_yticklabels(): tl.set_color(col1)
for tl in ax2.get_yticklabels(): tl.set_color(col2)

ax2.set_yscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')

labels = [l.get_label() for l in plots]
plt.legend(plots,labels,numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)

#plt.subplots_adjust(bottom=0.2,left=0.1,right=0.95)

# SAVE
plt.close()

