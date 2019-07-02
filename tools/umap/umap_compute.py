import sys
args = sys.argv

pca_data_fname = args[1]
umap_coordinates_fname = args[2]

import pandas as pd
import umap
    
print("Loading data ...")
pca = pd.read_csv(pca_data_fname, sep = ",", index_col = 0).values

print("creating UMAP ... ")
embedding = umap.UMAP(n_neighbors = 5, min_dist = .3, metric = "correlation", random_state = 17)
embedding = embedding.fit_transform(pca)

print("saving UMAP coordinates to disk ...")
df = pd.DataFrame(embedding)
df.columns = ["UMAPx", "UMAPy"]
df.to_csv(umap_coordinates_fname)
