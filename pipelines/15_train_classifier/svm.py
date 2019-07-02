import pandas as pd
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
import pickle

import sys
args = sys.argv
material_dir = args[1]
output_dir = args[2]

from os.path import join

print("Loading data ...")
X = pd.read_csv(join(material_dir, 'data.csv'), sep = ",", index_col = 0).values
Y = pd.read_csv(join(material_dir, 'labels.csv')).values[:, 0].reshape(-1, 1).ravel()

from sklearn.decomposition import PCA
pca = PCA(n_components = .8)
X = pca.fit_transform(X)

modelFile = open(join(output_dir, "pca.pickle"), "wb")
print(modelFile)
modelFile.write(pickle.dumps(pca))
modelFile.close()

print("Splitting into training and test sets...")
(X_train, X_test, y_train, y_test) = train_test_split(X, Y, test_size = .3, random_state = 42)

params = {"C":[1e-6, 1e-3, .1, 1, 10, 100, 1000],
          "gamma": [1e-6, 1e-3, .1, 1]}

# established as the best paramaters in some other work
params = {"C":[10], "gamma": [1e-3]}

print("Creating the model and fitting the data ...")
model = GridSearchCV(SVC(probability = False, kernel = "rbf"), params, cv=5)
model.fit(X_train, y_train)

print("Testing ...")
pred = model.predict(X_test)
cls_report = classification_report(y_test, pred, target_names = model.classes_)
print(cls_report)
with open(join(output_dir, 'classification_report.txt'), "w")  as cl_f:
    cl_f.write(cls_report)

print("Saving model and confusion matrix to disk ...")
cnf_matrix = confusion_matrix(y_test, pred)
df = pd.DataFrame(cnf_matrix)
df.columns = model.classes_
df.to_csv(join(output_dir, 'confusion_matrix.csv'))

modelFile = open(join(output_dir, 'model.pickle'), "wb")
modelFile.write(pickle.dumps(model))
modelFile.close()

