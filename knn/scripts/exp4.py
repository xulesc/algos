# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 12:12:31 2014

@author: anuj
"""

print(__doc__)

import numpy as np

from sklearn.preprocessing import normalize
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cross_validation import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn import preprocessing
from sklearn.naive_bayes import MultinomialNB

d_dir = '/home/xule/workspace_git/uni_courses/phd_courses/pms_510'
klas_file = '%s/assignment/data/diabetic_data.csv.klass' %d_dir
data_file = '%s/assignment/data/diabetic_data.csv.non_num' %d_dir
data_file2 = '%s/assignment/data/diabetic_data.csv.num' %d_dir

def metric_from_cm(cm):
    a = cm[1][1]; b = cm[0][1]; c = cm[1][0]; d = cm[0][0]
    return [1.0 * a / (a + b + .00001), 1.0 * a / (a + c + .00001), 1.0 * d / (d + b + .00001)]

def klassify(name, f_data, labels, grade = True, neigh = MultinomialNB()):
    print '#dataset size: %d' %len(f_data)
    if grade:
        data_train, data_test, labels_train, labels_test = train_test_split(f_data, labels, test_size=0.20, random_state=42)
    else:
        data_train = f_data
        labels_train = labels
    #neigh = KNeighborsClassifier(n_neighbors=1)
    #from sklearn.naive_bayes import MultinomialNB
    #neigh = MultinomialNB()
    neigh.fit(data_train, labels_train)
    
    if grade == False:
        return neigh
        
    pred = neigh.predict(data_test)        
    pre, rec, tnr = metric_from_cm(confusion_matrix(labels_test, pred))
    print('%s: %f (pre), %f (rec), %f (tnr)\n' %(name, pre, rec, tnr))
    
def klassify2(name, f_data_num, f_data_non_num, labels):    
    thresh = int(1.0 * 80 * len(labels) / 100)
    num_neigh = klassify('t', f_data_num[:thresh], labels[:thresh], False, neigh = KNeighborsClassifier(n_neighbors=1))
    non_num_neigh = klassify('t', f_data_non_num[:thresh], labels[:thresh], False)
    #data_p = zip(num_neigh.predict(f_data_num[:thresh]), non_num_neigh.predict(f_data_non_num[:thresh]))
    #neigh = KNeighborsClassifier(n_neighbors=1)
    #neigh.fit(data_p, labels[:thresh])
    #pred = neigh.predict(zip(num_neigh.predict(f_data_num[thresh:]), non_num_neigh.predict(f_data_non_num[thresh:])))        
    num_p = num_neigh.predict(f_data_num[thresh:])
    non_num_p = non_num_neigh.predict(f_data_non_num[thresh:])
    pred = [int(a) or int(b) for a,b in zip(num_p, non_num_p)]
    pre, rec, tnr = metric_from_cm(confusion_matrix(labels[thresh:], pred))
    print('%s: %f (pre), %f (rec), %f (tnr)\n' %(name, pre, rec, tnr))

##############
klasses = np.genfromtxt(klas_file, skip_header=1)
n_data = np.genfromtxt(data_file, skip_header=1,delimiter=',',dtype='|S5')
n_data_num = np.genfromtxt(data_file2, usecols=range(2, 11), skip_header=1,delimiter=',', missing_values='?')

exc = np.isnan(n_data_num).any(axis=1)
n_data_num_n = normalize(n_data_num[~exc], axis = 0)
labels = klasses[~exc]
n_data2 = n_data[~exc]
n_data2 = [x[:len(x) - 1] for x in n_data2]

n_data2 = np.transpose(n_data2)
le = preprocessing.LabelEncoder()
n_data3 = [le.fit(d).transform(d) for d in n_data2]
##############
f_data = np.transpose(n_data3)
klassify('non-numeric', f_data, labels)
f_data = n_data_num_n
klassify('numeric', f_data, labels, neigh = KNeighborsClassifier(n_neighbors=1))
f_data = n_data3
for x in np.transpose(n_data_num_n):
    f_data.append(x)
f_data = np.transpose(f_data)
klassify('mixed', f_data, labels)
##############
klassify2('ensemble', n_data_num_n, np.transpose(n_data3), labels)

##############


##

