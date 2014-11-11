# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 10:14:30 2014

@author: xule
"""
import csv

import numpy as np

from sklearn import preprocessing
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.lda import LDA

from sklearn.feature_selection import SelectFwe, SelectKBest, f_regression, chi2
from sklearn.pipeline import Pipeline

TRG_SIZE = 3200
##
csv_reader = lambda x: csv.reader(open(x, 'rb'), delimiter=',')
base_stats = lambda x: [min(x), max(x), np.mean(x), np.std(x)]
min_max_scaler = preprocessing.MinMaxScaler()
##

data_file='/home/xule/workspace/pms_510/assignment/data/abalone.data'
header=['sex','length','diameter','height','whole','shucked','viscera','shell','rings']
data=zip(*map(lambda x: map(lambda y: float(y), x), csv_reader(data_file)))
## Stats for numeric domains
stats = map(base_stats, data)
#print stats
##
names = ["Nearest Neighbors",  "Linear SVM", "RBF SVM", 
         "Decision Tree", "Random Forest", "AdaBoost", 
         "Naive Bayes", "LDA"] #, "QDA"]         
anova_filter = SelectFwe(chi2) #SelectKBest(chi2, k=4)
classifiers = [ KNeighborsClassifier(3), SVC(kernel="linear", C=0.025),
    SVC(gamma=2, C=1), DecisionTreeClassifier(max_depth=5),
    RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
    AdaBoostClassifier(), GaussianNB(), LDA()] #, QDA()]
##
X = zip(*data[:len(data)-1]); Y = data[len(data)-1]
X_trg = min_max_scaler.fit_transform(X[:TRG_SIZE]); Y_trg = Y[:TRG_SIZE]
X_tst = min_max_scaler.fit_transform(X[TRG_SIZE:]); Y_tst = Y[TRG_SIZE:]
for name, classifier in zip(names, classifiers):
    classifier.fit(X_trg, Y_trg)
    print '%s: %f' %(name, classifier.score(X_tst, Y_tst))
    clf = Pipeline([('anova', anova_filter), (name, classifier)])
    clf.fit(X_trg, Y_trg)
    print '%s: %f' %('anova-%s' %name, clf.score(X_tst, Y_tst))    
    
##
