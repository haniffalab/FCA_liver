#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 15:01:36 2018

@author: doru
"""

import sys
args = sys.argv
working_folder = args[1]

import matplotlib; matplotlib.use('Agg');
import scanpy.api as sc;
import pandas as pd

sc.settings.verbosity = 3

scObj = sc.read("{CW}/raw_data.mtx".format(CW = working_folder), cache = False).T

# load gene names
scObj.var_names = pd.read_csv("{CW}/genenames.csv".format(CW = working_folder)).iloc[:, 1]

# load cell names
scObj.obs_names = pd.read_csv("{CW}/cellnames.csv".format(CW = working_folder)).iloc[:, 1]

# add cell labels
cell_labels = pd.read_csv("{CW}/cell_labels.csv".format(CW = working_folder), index_col = 0)

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

# save diffusion map to disk
dm  = scObj.obsm["X_diffmap"]
dm = pd.DataFrame(data = dm, index = None, columns = None)
dm.to_csv("{CW}/dm.csv".format(CW = working_folder), columns = None, header = None)


