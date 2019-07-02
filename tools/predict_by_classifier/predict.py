import sys
args = sys.argv
model_fname = args[1]
data_fname = args[2]
predictions_fname = args[3]
pca_fname = args[4]

import pickle
import pandas as pd

print("loading model ...")
model = pickle.loads(open(model_fname, "rb").read())
pca   = pickle.loads(open(pca_fname, "rb").read())

print("loading data....")
predictions = pd.read_csv(data_fname, sep = ",", index_col = 0).values
print("classifying cells ... ")
predictions = pca.transform(predictions);
predictions = model.predict(predictions)
print("writing results to disk...")
df = pd.DataFrame(predictions)
df.to_csv(predictions_fname)
print("Classification finished. Python will now say goodbye. Handing over to R.")