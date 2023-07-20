#!/usr/bin/env python3
#
#   https://scikit-learn.org/stable/modules/outlier_detection.html#outlier-detection
#
#   cat color-x1-x2-x3-....txt   | mplot reduce -o img.png
#
#   color .. if RGB string, escape the hash as \#00ff00
#
# CMDLINE

from collections import OrderedDict
from functools import partial
from time import time
from sklearn import manifold

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

files  = []
# FILES
fname = files[0]

# # Next line to silence pyflakes. This import is needed.
# Axes3D

colors = []
dat = []
reader = open(fname)
for line in reader:
    row = line.rstrip('\n').split('\t')
    if row[0][0]=='\\': row[0] = row[0][1:]
    colors.append(row[0])
    row = row[1:]
    for i in range(len(row)): row[i] = float(row[i])
    dat.append(row)

#   n_points = 100
#   X, color = datasets.make_s_curve(n_points, random_state=0)
n_neighbors = 10
n_components = 2

# Create figure
fig = plt.figure(figsize=(10, 8))

# Set-up manifold methods
LLE = partial(manifold.LocallyLinearEmbedding, n_neighbors=n_neighbors, n_components=n_components, eigen_solver='auto')

methods = OrderedDict()
methods['LLE'] = LLE(method='standard')
methods['LTSA'] = LLE(method='ltsa')
methods['Hessian LLE'] = LLE(method='hessian')
#methods['Modified LLE'] = LLE(method='modified')
#methods['Isomap'] = manifold.Isomap(n_neighbors, n_components)
#methods['MDS'] = manifold.MDS(n_components, max_iter=100, n_init=1)
#methods['SE'] = manifold.SpectralEmbedding(n_components=n_components, n_neighbors=n_neighbors)
#methods['t-SNE'] = manifold.TSNE(n_components=n_components, init='pca', random_state=0)

# Plot results
for i, (label, method) in enumerate(methods.items()):
    t0 = time()
    Y = method.fit_transform(dat)
    t1 = time()
    print("%s: %.2g sec" % (label, t1 - t0))
    ax = fig.add_subplot(2, 5, 2 + i + (i > 3))
    ax.scatter(Y[:, 0], Y[:, 1], c=colors)
    ax.set_title("%s (%.2g sec)" % (label, t1 - t0))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.axis('tight')
    break;

# dpi: dpi

# SAVE
plt.close()


