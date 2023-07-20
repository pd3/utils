#!/usr/bin/env python3
#
# mplot gene -o rmme.png +gff $gff +g MEF2C variants.txt
#
#   variants.txt
#       DEL chr beg end
#       SNP chr pos color
#
#   options
#       +line 'pos;color=red;lw=2;ls=:
#       +track file,color=red       # file: chr,beg,end,value;      e.g. BroadH3K4Me3_MEF2Cregion.txt
#       +cl utr=red                 # override default colors
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys,subprocess,re

def error(msg):
    print(msg)
    sys.exit(1)

files  = []
# FILES

gff = None
# gff: 'gff'
if gff==None: error('Missing the +gff option')

gene = None
# g: 'gene'
if gene==None: error('Missing the +g option')

def run_cmd(cmd):
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return proc.stdout

def parse_gff(gff,gene):
    dat = {}
    gene_id = None
    tscript_id = None
    tscripts = {}
    if re.search(r'.gz$',gff):
        cmd = "gunzip -c "+gff
    else:
        cmd = cmd = 'cat '+gff
    for line in run_cmd(cmd):
        line = line.decode('ascii')
        if line[0]=='#': continue
        line.rstrip('\n')
        col = line.split('\t')
        if col[2]=="gene":
            if gene_id!=None: break
            if col[8].find('Name='+gene+';')==-1: continue
            m = re.search('ID=gene:([^;]+)',col[8])
            gene_id = m.group(1)
            dat['chr'] = col[0]
            dat['beg'] = float(col[3])/1e6
            dat['end'] = float(col[4])/1e6
            if dat['beg'] > dat['end']: tmp = dat['beg']; dat['beg'] = dat['end']; dat['end'] = tmp;
            dat['strand'] = col[6]
        if gene_id==None: continue
        if col[8].find('Parent=gene:'+gene_id+';')!=-1:
            if col[8].find('biotype=protein_coding')==-1:
                tscript_id = None
                continue
            m = re.search('ID=transcript:([^;\s]+)',col[8])
            if m==None: continue
            tscript_id = m.group(1)
        if tscript_id==None: continue
        m = re.search('Parent=transcript:([^;\s]+)',col[8])
        if m==None: continue
        if col[2]=="CDS":
            tscript_id = m.group(1)
            if tscript_id not in tscripts: tscripts[tscript_id] = []
            beg = float(col[3])/1e6
            end = float(col[4])/1e6
            if beg > end: tmp = beg; beg = end; end = tmp
            tscripts[tscript_id].append([beg,end,'CDS'])
        if col[2]=="three_prime_UTR" or col[2]=="five_prime_UTR":
            tscript_id = m.group(1)
            if tscript_id not in tscripts: tscripts[tscript_id] = []
            beg = float(col[3])/1e6
            end = float(col[4])/1e6
            if beg > end: tmp = beg; beg = end; end = tmp
            tscripts[tscript_id].append([beg,end,'UTR'])
    dat['tscripts'] = tscripts
    return dat

def parse_variants(fname):
    vars = {}
    f = open(fname)
    for line in f:
        line = line.rstrip("\n")
        col = line.split('\t')
        if col[0]=='DEL' or col[0]=='DUP':
            type = col[0]
            if type not in vars: vars[type] = []
            beg = float(col[2])/1e6
            end = float(col[3])/1e6
            if beg > end: tmp = beg; beg = end; end = tmp
            vars[type].append([beg,end])
        elif col[0]=='SNP':
            if 'SNP' not in vars: vars['SNP'] = []
            color = '#D43F3A'
            if len(col)==4: color = col[3]
            vars['SNP'].append([float(col[2])/1e6,color])
    return vars

vars = parse_variants(files[0])
dat  = parse_gff(gff,gene)
tscripts = dat['tscripts']

wh = (10,5)
# wh: wh
fig, ax1 = plt.subplots(1, 1, figsize=wh)

ht = 0.3
# ht: ht

color_tr  = '#242424'
color_utr = '#f4640d'
cl = None               # +cl utr=red
# cl: 'cl'
if cl!=None:
    tmp = cl.split(',')
    for x in tmp:
        key,col = x.split('=')
        if key=='utr': color_utr = col
        else: error('todo: color for '+key)

dtrack = 0.3
dtscript  = 0.5
itscript  = 0
for tr in tscripts:
    itscript += dtscript
    tscript = tscripts[tr]
    for i in range(len(tscript)-1):
        ax1.plot([tscript[i][1],tscript[i+1][0]],[itscript+ht*0.5,itscript+ht*0.5],color=color_tr,zorder=0)
    for ex in tscript:
        color = 'black'
        if ex[2]=='UTR': color = color_utr
        rect = patches.Rectangle((ex[0],itscript), ex[1]-ex[0], ht, color=color,zorder=10)
        ax1.add_patch(rect)

delta = (dat['end']-dat['beg'])*0.1
xmin = dat['beg'] - delta
xmax = dat['end'] + delta

# Add arrows to show the strand
itscript  = 0
for tr in tscripts:
    itscript += dtscript
    tscript = tscripts[tr]
    imax = 0 
    for i in range(len(tscript)-1):
        if tscript[imax+1][0]-tscript[imax][1] < tscript[i+1][0]-tscript[i][1]: imax = i
    if len(tscript)>1:
        mid  = (tscript[imax+1][0]+tscript[imax][1])*0.5
        alen = (xmax-xmin)*0.01
        if dat['strand']=='-': alen = -alen
        ax1.annotate("",xy=(mid,itscript+ht*0.5),xytext=(mid-alen,itscript+ht*0.5),arrowprops=dict(color=color_tr,headwidth=7),zorder=0)


# custom tracks
usr = None
# track: ['usr']
if usr!=None:
    for track in usr:
        tmp = track.split(',')
        tmp_args = {}
        for x in tmp[1:]:
            a,b = x.split('=')
            tmp_args[a] = b
        f = open(tmp[0])
        x = []
        y = []
        max = None
        for line in f:
            if line[0]=='#': continue
            line = line.rstrip("\n")    # chr,beg,end,value
            val  = line.split('\t')
            x.append(float(val[1])/1e6)
            y.append(float(val[3]))
            if max==None or max<y[-1]: max = y[-1]
        itscript += dtscript
        ax1.fill_between(x,[ht*v/max+itscript for v in y],itscript,**tmp_args)

# SNPs
itrack = 0
if 'SNP' in vars:
    import random
    xdelta = (xmax-xmin)/100
    ytrack = {}
    for x in  vars['SNP']:
        xapprox = xmin + int((x[0]-xmin)/xdelta) * xdelta
        if xapprox not in ytrack: ytrack[xapprox] = 0
        y = ytrack[xapprox]
        if itrack < y: itrack = y
        ytrack[xapprox] += 1
        yjit = y + random.random()*0.6 - 0.3
        ax1.plot([x[0]],[-yjit],'*',color=x[1])

# Deletions
if 'DEL' in vars:
    for i in range(len(vars['DEL'])):
        delv = vars['DEL'][i]
        if delv[1] < xmin or delv[0] > xmax: continue
        rect = patches.Rectangle((delv[0],-dtrack*itrack), delv[1]-delv[0], dtrack*0.5, color='#faa300',zorder=10)
        ax1.add_patch(rect)
        itrack += 1

# Duplications
if 'DUP' in vars:
    for i in range(len(vars['DUP'])):
        dupv = vars['DUP'][i]
        if dupv[1] < xmin or dupv[0] > xmax: continue
        rect = patches.Rectangle((dupv[0],-dtrack*itrack), dupv[1]-dupv[0], ht, color='#4CAE4C',zorder=10)
        ax1.add_patch(rect)
        itrack += 1

lines = None         # +line 'pos;color=red;lw=2;ls=:;marker=o;ms=10'  (can be given multiple times)
# line: ['lines']
if lines!=None:
    for line in lines:
        tmp = line.split(';')
        lin = tmp[0].split(',')
        tmp_args = {}
        for x in tmp[1:]:
            a,b = x.split('=')
            tmp_args[a] = b
        ax1.plot([float(lin[0])/1e6,float(lin[0])/1e6],[0.5,-dtrack*(itrack-1)],**tmp_args)

title = None
if dat['strand']=='-': title = '< '+gene
else: title = gene+' >'

ax1.set_xlim(xmin,xmax)
ax1.set_ylim(-itrack*dtrack,itscript+dtscript)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.get_xaxis().tick_bottom()
ax1.set_yticks([])
ax1.spines['bottom'].set_color('grey')
ax1.patch.set_visible(False)
ax1.set_title(title)
ax1.set_xlabel('chr'+dat['chr']+' [Mb]')
ax1.xaxis.labelpad = 10

adjust = None
# adj: {adjust}
if adjust!=None: plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95

# dpi: dpi
# SAVE
plt.close()

