import pandas as pd
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
import pickle

import sys
from os.path import join
args = sys.argv
save_to = args[1]

print("Loading data ...")
X = pd.read_csv(join(save_to, "./data.csv"), sep = ",", index_col = 0).values
y = pd.read_csv(join(save_to, "./labels.csv")).values[:, 0].reshape(-1, 1).ravel()

from sklearn.decomposition import PCA
pca = PCA(n_components = .8)
X = pca.fit_transform(X)

modelFile = open(join(save_to, "pca.pickle"), "wb")
print(modelFile)
modelFile.write(pickle.dumps(pca))
modelFile.close()

print(X.shape)

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=19)

params = {"C":[1, 10, 100, 300]}

print("Creating the model and fitting the data ...")
model = GridSearchCV(SVC(probability = False, kernel = "rbf"), params, cv=3)
model.fit(X_train, y_train)

print("Testing ...")
pred = model.predict(X)
cls_report = classification_report(y, pred, target_names = model.classes_)
print(cls_report)
with open(join(save_to, "classification_report.txt"), "w")  as cl_f:
    cl_f.write(cls_report)

print("Saving model and confusion matrix to disk ...")
cnf_matrix = confusion_matrix(y, pred)
df = pd.DataFrame(cnf_matrix)
df.columns = model.classes_
df.to_csv(join(save_to, "confusion_matrix.csv"))

modelFile = open(join(save_to, "model.pickle"), "wb")
modelFile.write(pickle.dumps(model))
modelFile.close()

print(model.best_params_)
