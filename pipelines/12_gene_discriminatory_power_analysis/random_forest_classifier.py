#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 11:17:19 2018

@author: doru
"""

import sys
from os.path import join

arguments = sys.argv

working_folder = arguments[1]
material_folder = join(working_folder, "material")

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.manifold import TSNE
import pandas as pd
import numpy as np
from numpy.core.umath_tests import inner1d

class randomGuessing:
    def __init__(self, Y_Labels):
        tally = np.unique(Y_Labels, return_counts = True)
        self.labels = list(tally[0])
        frequency = tally[1]
        frequency = frequency / frequency.sum()
        self.frequency = frequency
    def predict(self, X):
        return np.array(np.random.choice(self.labels, X.shape[0], p = self.frequency))

print("Loading data ...")
X = pd.read_csv(join(material_folder, "data.csv"), sep = ",", index_col = 0).values
labels_and_colours = pd.read_csv(join(material_folder, "labels.csv")).values
Y = labels_and_colours[:, 0].reshape(-1, 1).ravel()
Colours = labels_and_colours[:, 1].reshape(-1, 1).ravel()
# separate labels and colours

print("Splitting into training and test sets...")
(X_train, X_test, y_train, y_test) = train_test_split(X, Y, test_size = .3, random_state = 32)

randomForestClassifier = RandomForestClassifier(n_estimators = 500, criterion = "gini", min_samples_split = 5, bootstrap = True, n_jobs=20,class_weight="balanced")

randomForestClassifier.fit(X_train, y_train)
pred = randomForestClassifier.predict(X_test)

cls_report = classification_report(y_test, pred, target_names = randomForestClassifier.classes_)
with open(join(working_folder, "classification_report.txt"), "w")  as cl_f:
    cl_f.write(cls_report)

print("Saving confusion matrix to disk ...")
cnf_matrix = confusion_matrix(y_test, pred)
df = pd.DataFrame(cnf_matrix)
df.columns = randomForestClassifier.classes_
df.to_csv(join(material_folder, "confusion_matrix.csv"))

print("Finishing random_forest_classifier.py run")
