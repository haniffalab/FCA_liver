import sys
args = sys.argv
n_iterations = int(args[1])
pca_data_fname = args[2]
snn_fname = args[3]
fdg_coordinates_fname = args[4]

from fa2 import ForceAtlas2
import pandas as pd
from scipy.io import mmread
import numpy as np

# load pca, SNN and label colours data
# the first 2 PC form PCA are used as initial conditions
# SNN is used for building the force directed graph
pca_data = pd.read_csv(pca_data_fname, index_col = 0)
snn      = mmread(snn_fname)

# set initialposition as the first 2 PCs
positions = pca_data.values[:, 0:2]

# initialize force directed graph class instance
forceatlas2 = ForceAtlas2(outboundAttractionDistribution=False, linLogMode=False,
                          adjustSizes=False, edgeWeightInfluence=1.0,
                          jitterTolerance=1.0, barnesHutTheta = .8,
                          barnesHutOptimize=True, multiThreaded=False,
                          scalingRatio=2.0, strongGravityMode=True, gravity=1, verbose=True)

# run force directed graph; for each iterations generates the coordinates use din each frame
coord = forceatlas2.forceatlas2(G = snn, pos = positions, iterations = n_iterations)
coord = np.array(coord)

np.savetxt( fname = fdg_coordinates_fname, X = coord, delimiter = ",")