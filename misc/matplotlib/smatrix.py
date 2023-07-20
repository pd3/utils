#!/usr/bin/python3.4
#
# similarity matrix
#
#   cat id-id-dist.txt | mplot smatrix -o smatrix.png
#

# CMDLINE

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import csv, random

from scipy.spatial.distance import pdist, squareform
from sklearn import datasets
from fastcluster import linkage

csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

def seriation(Z,N,cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z
            
        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index-N,0])
        right = int(Z[cur_index-N,1])
        return (seriation(Z,N,left) + seriation(Z,N,right))
    
def compute_serial_matrix(dist_mat,method="ward"):
    '''
        input:
            - dist_mat is a distance matrix
            - method = ["ward","single","average","complete"]
        output:
            - seriated_dist is the input dist_mat,
              but with re-ordered rows and columns
              according to the seriation, i.e. the
              order implied by the hierarchical tree
            - res_order is the order implied by
              the hierarhical tree
            - res_linkage is the hierarhical tree (dendrogram)
        
        compute_serial_matrix transforms a distance matrix into 
        a sorted distance matrix according to the order implied 
        by the hierarchical tree (dendrogram)
    '''
    N = len(dist_mat)
    flat_dist_mat = squareform(dist_mat)
    res_linkage = linkage(flat_dist_mat, method=method,preserve_input=True)
    res_order = seriation(res_linkage, N, N + N-2)
    seriated_dist = np.zeros((N,N))
    a,b = np.triu_indices(N,k=1)
    seriated_dist[a,b] = dist_mat[ [res_order[i] for i in a], [res_order[j] for j in b]]
    seriated_dist[b,a] = seriated_dist[a,b]
    return seriated_dist, res_order, res_linkage


type = 'iid'
# type: 'iid'

files = []
# FILES
file = files[0]

(xmin,xmax) = (None,None)
xlim = None
# xr: xlim
if xlim!=None:
    (xmin,xmax) = xlim

mat = []
raw = {}
ids = {}
with open(file, 'r') as f:
    reader = csv.reader(f, 'tab')
    for row in reader:
        if row[0][0] == '#': continue
        id1  = row[0]
        id2  = row[1]
        dist = row[2]
        if id1 not in ids: ids[id1] = len(ids)
        if id2 not in ids: ids[id2] = len(ids)
        i1 = ids[id1]
        i2 = ids[id2]
        key = "%d.%d" % (i1,i2)
        raw[key] = [i1,i2,float(dist)]
    mat = np.zeros(shape=(len(ids),len(ids)))
    for x in raw.values():
        mat[x[0]][x[1]] = x[2]
        mat[x[1]][x[0]] = x[2]
N = len(ids)

# np.random.seed(seed=1)
# iris = datasets.load_iris()
# N = len(iris.data)
# X = iris.data[np.random.permutation(N),:]
# mat = squareform(pdist(X))

method = 'average'     # one of: ward single average complete centroid median;  e.g. +mt ward
# mt: 'method'
ordered_dist_mat, res_order, res_linkage = compute_serial_matrix(mat,method)
    
plt.pcolormesh(ordered_dist_mat)
plt.xlim([0,N])
plt.ylim([0,N])

# dpi: dpi

# SAVE
plt.close()

