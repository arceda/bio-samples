# this script compare the performace of:
# Toward an Alignment-Free Method for Feature Extraction [1]
# An open-source k-mer based machine learning tool for fast and accurate [2]

import matplotlib.pyplot as plt
from sklearn.metrics import f1_score
from sklearn.feature_selection import RFE
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_recall_fscore_support

from sklearn.svm import SVC
from sklearn.naive_bayes import MultinomialNB
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
import time
import feature_extractor as fe
from random import shuffle
import sys
import numpy as np

import mykameris as kam

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

###################################################################################################
#############################               kameris                 ###############################
# this funtion return de acc, fcore, .... of the models proposed by
# An open-source k-mer based machine learning tool for fast and accurate [2]
def kameris(train, test, k):
    X_train = []
    y_train = []
    for sample in train:
        y_train.append(sample[2])
        seq = sample[1]
        k_mers_frecuencies = kam.cgr(seq, k)
        k_mers_frecuencies = k_mers_frecuencies/sum(k_mers_frecuencies) #normalize
        X_train.append(k_mers_frecuencies)  
    X_train = np.matrix(X_train)
    #print(X_train.shape, X_train)
    #print(y_train)

    clf = SVC(kernel = "linear", C = 1) 
    clf.fit(X_train, y_train)

    X_test = []
    y_test = []
    for sample in test:
        y_test.append(sample[2])
        seq = sample[1]
        k_mers_frecuencies = kam.cgr(seq, k)
        k_mers_frecuencies = k_mers_frecuencies/sum(k_mers_frecuencies) #normalize
        X_test.append(k_mers_frecuencies)
    X_test = np.matrix(X_test)

    y_pred = clf.predict(X_test)
        
    metrics = precision_recall_fscore_support(y_test, y_pred, average='weighted')
    acc = accuracy_score(y_test, y_pred)

    #print('metrics: acc, precision, recall, fscore ', acc, metrics)
    return acc, metrics[0], metrics[1], metrics[2]


###################################################################################################
#############################               castor                  ###############################
# this funtion return de acc, fcore, .... of the models proposed by
# Toward an Alignment-Free Method for Feature Extraction [1]
def castor(train, test, k, nfeatures, rfe = False):
    k_mers              = fe.generate_K_mers(train, k) #list of substring of size k: (if k = 2; k_mers= [AT, CG, AC, ...])    
    _k_mers             = k_mers.copy()
    X_train, y_train    = fe.generateXYMatrice(train, k_mers, k) # OCURRENCE MATRIX
    X_train             = fe.maxMinNormalization(X_train)

    if rfe:
        X_train, k_mers   = fe.recursiveFeatureElimination(X_train, y_train, k_mers, nfeatures)
    
    # Implement and fit classifier
    clf = SVC(kernel = "linear", C = 1) 
    clf.fit(X_train, y_train)

    # test
    X_test, y_test      = fe.generateXYMatrice(test, k_mers, k) # OCURRENCE MATRIX
    X_test              = fe.maxMinNormalization(X_test)

    y_pred = clf.predict(X_test)
        
    metrics = precision_recall_fscore_support(y_test, y_pred, average='weighted')
    acc = accuracy_score(y_test, y_pred)

    #print('metrics: acc, precision, recall, fscore ', acc, metrics)
    return acc, metrics[0], metrics[1], metrics[2]

    

###################################################################################################
###################################################################################################

def split_data(data, indexs):
    k_fold = []
    for i in indexs:
        k_fold.append(data[i])
    
    return k_fold

#k = 7
nfeatures = 7

datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 
            'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT'] 

dataset_path = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/"
#dataset_path = sys.argv[1]

clf = SVC(kernel = "linear", C = 1)  


for i, dataset in enumerate(datasets):

    print(i, "\n\nEVALUATING DATASET: ", dataset, "...")
    print('===========================================================')
    print('===========================================================') 

    for k in range(1,11):
        print("\n\nEvaluating with k-mer:", k, " ==========================")

        data = fe.generateLabeledData(dataset_path + dataset + "/data.fa", dataset_path  + dataset + "/class.csv")         
        
        kf = KFold(n_splits=5, shuffle=True, random_state=1)
        i = 0

        metrics_castor = []
        metrics_kameris = []
        for train_index, test_index in kf.split(np.zeros(len(data))):            

            data_train = split_data(data, train_index)
            data_test = split_data(data, test_index)                       

            acc, pre, recall, fscore = kameris(data_train, data_test, k)
            metrics_kameris.append([acc, pre, recall, fscore])
            print("k-fold ", i, "metrics_kameris: ", acc, pre, recall, fscore)

            acc, pre, recall, fscore = castor(data_train, data_test, k, nfeatures)
            metrics_castor.append([acc, pre, recall, fscore])
            print("k-fold ", i, "metrics_castor: ", acc, pre, recall, fscore)

            i += 1
            
        metrics_kameris = np.matrix(metrics_kameris)
        metrics_castor = np.matrix(metrics_castor)
        print("mean metrics_kameris: ", metrics_kameris.mean(0))
        print("mean metrics_castor: ", metrics_castor.mean(0))
        
    break

