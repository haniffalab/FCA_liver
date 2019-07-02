#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 03 10:49:22 2018

@author: doru
"""

import sys
args = sys.argv
output_folder = args[1]
n_clusters = int(args[2])

import os

os.chdir(output_folder)

from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph
from sklearn import mixture
from sklearn import metrics
import pandas as pd
import numpy as np

print("Loading data ...")
pca_df = pd.read_csv("material/pca.csv", sep = ",", index_col = 0)
pca = pca_df.values[:, :20]
ori_labels = pca_df.values[:, 20]

print("computing agglomerative clustering ... ")
connectivity_graph = kneighbors_graph(X = pca, n_neighbors = 30)
clustering = AgglomerativeClustering(n_clusters = n_clusters, affinity = "euclidean", 
                                     connectivity = connectivity_graph, 
                                     linkage = "ward")
agg_labels = clustering.fit(pca)
agg_labels = agg_labels.labels_

print("computing gaussian mixture clustering ... ")
gmm = mixture.GaussianMixture(n_components = n_clusters, covariance_type = "full",
                              random_state = 52).fit(pca)
gmm_labels = gmm.predict(pca)

# metrics

# voted labels
def voted_labels(true, pred):
    pred_unique = np.unique(pred)
    pred_voted = true.copy()
    for p in pred_unique:
        p_indx = pred == p
        ori_types = true[p_indx]
        ori_unique = np.unique(ori_types, return_counts = True)
        voted = ori_unique[0][np.argmax(ori_unique[1])]
        pred_voted[p_indx] = voted
    return pred_voted

agg_labels = voted_labels(ori_labels, agg_labels)
gmm_labels = voted_labels(ori_labels, gmm_labels)

# compute adjusted Rand Index
# reference: Comparing partitions - https://link.springer.com/article/10.1007%2FBF01908075
RandIndex_Louvain_Agg = metrics.adjusted_rand_score(ori_labels, agg_labels)
RandIndex_Louvain_GMM = metrics.adjusted_rand_score(ori_labels, gmm_labels)

# ref for mutual information metrics - "Information theoretic measures for clusterings comparison"

# compute adjusted mutual information
AdjMutInf_Louvain_Agg = metrics.adjusted_mutual_info_score(ori_labels, agg_labels)
AdjMutInf_Louvain_GMM = metrics.adjusted_mutual_info_score(ori_labels, gmm_labels)

RandIndex_Louvain_Agg = str(RandIndex_Louvain_Agg)
RandIndex_Louvain_GMM = str(RandIndex_Louvain_GMM)
AdjMutInf_Louvain_Agg = str(AdjMutInf_Louvain_Agg)
AdjMutInf_Louvain_GMM = str(AdjMutInf_Louvain_GMM)

# save measures to disk
result = '\n'.join([RandIndex_Louvain_Agg, RandIndex_Louvain_GMM, 
                    AdjMutInf_Louvain_Agg, AdjMutInf_Louvain_GMM])

with open('material/agreement_measures.txt', 'w') as agg_file:
  agg_file.write(result)
  
print("saving gaussian mixture clustering labels to disk ...")
df = pd.DataFrame(gmm_labels)
df.index = pca_df.index
df.to_csv("./material/gaussian_mixture.csv")

print("saving agglomerative clustering labels to disk ...")
df = pd.DataFrame(agg_labels)
df.index = pca_df.index
df.to_csv("./material/agglomerative_clustering.csv")
