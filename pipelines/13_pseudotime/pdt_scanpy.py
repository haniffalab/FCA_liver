#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 15:01:36 2018

@author: doru
"""
print("starting .py script")

import sys
args = sys.argv
root_cell_type = args[1]
CWD = args[2]
print("printing args")
print(args)
args
# use the args below if you have a root cell type containing spaces and @@'s
#root_cell_type = args[1] + " " + args[2]
#CWD = args[3]

import matplotlib; matplotlib.use('Agg');
import scanpy.api as sc;
import pandas as pd
import numpy as np

print("printing root_cell_type")
print(root_cell_type)
print("printing CWD")
print(CWD)

sc.settings.verbosity = 3

scObj = sc.read("{CWD}/material/raw_data.mtx".format(CWD=CWD), cache = False).T

# load gene names
scObj.var_names = pd.read_csv("{CWD}/material/genenames.csv".format(CWD=CWD)).iloc[:, 1]

# load cell names
scObj.obs_names = pd.read_csv("{CWD}/material/cellnames.csv".format(CWD=CWD)).iloc[:, 1]

# add cell labels
cell_labels = pd.read_csv("{CWD}/material/cell_labels.csv".format(CWD=CWD), index_col = 0)

scObj.obs["cell_labels"] = cell_labels

# filter out genes present in less than 3 cells
sc.pp.filter_genes(scObj, min_cells=3)

# log-normalize the data
scObj.raw = sc.pp.log1p(scObj, copy=True)
sc.pp.normalize_per_cell(scObj, counts_per_cell_after=1e4)

# variable genes
filter_result = sc.pp.filter_genes_dispersion(
        scObj.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
# subset data on variable genes
scObj = scObj[:, filter_result.gene_subset]
# not sure?
sc.pp.log1p(scObj)

# scale the data
sc.pp.scale(scObj, max_value=10)

# run pca
sc.tl.pca(scObj)

# compunte neighborhood graph
sc.pp.neighbors(scObj, n_neighbors = 15, n_pcs = 20, knn = True, random_state = 10, method = "gauss")

# compute diffusion map
sc.tl.diffmap(scObj, n_comps = 20)

# set root
scObj.uns['iroot'] = np.flatnonzero(scObj.obs['cell_labels'] == root_cell_type)[0]

# compute dpt
print("computing sc.tl.dpt")
sc.tl.dpt(scObj, n_dcs = 20)

# pdt is at scObj.obs["dpt_pseudotime"]
print("displaying pdt table stored in scObj")
print(scObj.obs["dpt_pseudotime"])
pdt = scObj.obs["dpt_pseudotime"].to_csv("{CWD}/material/pseudotime.csv".format(CWD=CWD))

# save the pseudotime
dm  = scObj.obsm["X_diffmap"]
dm = pd.DataFrame(data = dm, index = None, columns = None)
dm.to_csv("{CWD}/material/dm.csv".format(CWD=CWD), columns = None, header = None)


