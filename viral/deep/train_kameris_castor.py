# this script train a CNN
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

from numpy.random import seed
seed(1)
import tensorflow as tf
tf.random.set_seed(1)


import pickle
from sklearn.model_selection import KFold 
from sklearn.model_selection import train_test_split
import collections
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import pylab
import math
import os
import sys
import cv2
import numpy as np

from Bio import SeqIO
import re
from cgr import save_fcgr
import glob

import tensorflow as tf

from tensorflow.keras import datasets, layers, models
from keras.models import Sequential
from keras.layers import Conv2D, MaxPooling2D, Dense, Flatten
from keras.utils import to_categorical
from keras.utils.vis_utils import plot_model
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import precision_recall_fscore_support
from keras.layers.normalization import BatchNormalization
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import confusion_matrix
from sklearn import svm
from sklearn.metrics import accuracy_score
from numpy import argmax
import matplotlib.pyplot as plt

import statistics
import pywt
from scipy.fftpack import fft
import pandas as pd
import csv

import mykameris as kam
import feature_extractor as fe

import time

def kameris(X_train, y_train, X_test, y_test, k, dimention_reduction, database_name):

    minMaxScaler = MinMaxScaler(feature_range=(0, 1), copy = False)
    scaler1 = StandardScaler()
    
    #############################################
    # compute k-mer frecuences 
    X_train_new = []
    for seq in X_train:       
        k_mers_frecuencies = kam.cgr(seq, k)   
        X_train_new.append(k_mers_frecuencies)        
    X_train = np.matrix(X_train_new)

    X_test_new = []
    for seq in X_test:       
        k_mers_frecuencies = kam.cgr(seq, k)   
        X_test_new.append(k_mers_frecuencies)        
    X_test = np.matrix(X_test_new)
    #############################################
    # scaling     
    #X_train = minMaxScaler.fit_transform(X_train)  # this lose variance   

    scaler1.fit(X_train)
    X_train = scaler1.transform(X_train)    
    
    #############################################
    # dimention reduction   
    number_features = int((np.count_nonzero(X_train > 0)*0.1) /len(X_train) )
    if dimention_reduction == 1 and pow(4, k) > number_features and number_features > 4:   

        #SVD
        svd = TruncatedSVD(n_components=number_features)
        rows, cols = X_train.shape
        svd.fit(X_train)
        X_train = svd.transform(X_train)
        print("SVD aplied ... X_train.shape: ", [rows, cols], X_train.shape)

    #############################################
    # train  
    clf = svm.SVC(kernel = "linear", C = 1) 
    clf.fit(X_train, y_train)

    pickle.dump(clf, open(current_dir + "/models/" + database_name + '-kameris.joblib', 'wb'))

    #############################################
    # scaling 
    #X_test = minMaxScaler.fit_transform(X_test) 
    X_test             = scaler1.transform(X_test)

    #############################################
    # dimention reduction 
    if dimention_reduction == 1 and pow(4, k) > number_features and number_features > 4:
        #X_test = pca1.transform(X_test)
        X_test = svd.transform(X_test)

    #############################################
    # predict
    y_pred = clf.predict(X_test)
        
    metrics = precision_recall_fscore_support(y_test, y_pred, average='weighted')
    acc = accuracy_score(y_test, y_pred)

    #print('metrics: acc, precision, recall, fscore ', acc, metrics)
    #return acc
    return acc, metrics[0], metrics[1], metrics[2]


###################################################################################################
#############################               castor                  ###############################
# this funtion return de acc, fcore, .... of the models proposed by
# Toward an Alignment-Free Method for Feature Extraction [1]
def castor(X_train, y_train, X_test, y_test, k, dimention_reduction, database_name):
    scaler1 = StandardScaler()

    train = []
    for i in range(X_train.shape[0]):
        train.append([0, X_train[i], y_train[i]])
    test = []
    for i in range(X_test.shape[0]):
        test.append([0, X_test[i], y_test[i]])


    k_mers              = fe.generate_K_mers(train, k) #list of substring of size k: (if k = 2; k_mers= [AT, CG, AC, ...])    
    _k_mers             = k_mers.copy()
    X_train, y_train    = fe.generateXYMatrice(train, k_mers, k) # OCURRENCE MATRIX    
    #X_train             = fe.maxMinNormalization(X_train) #use standar scaler instea of maxminnormalization    

    #############################################
    # dimention reduction
    if dimention_reduction == 1:
        #print(k_mers)

        #X_train, k_mers   = fe.recursiveFeatureElimination(X_train, y_train, k_mers, nfeatures)
        len_kamers_before_rfe = len(k_mers)
        k_mers, best_k_length = fe.getBestKmersAndFeaturesMini(train, k, range(1,100), 0.99)
        print("len_kamers_before_rfe ", len_kamers_before_rfe, " reduce fetures ", len(k_mers), k_mers)
        X_train, y_train    = fe.generateXYMatrice(train, k_mers, k)

   
        
    
    #############################################
    # scaling
    scaler1.fit(X_train)
    X_train             = scaler1.transform(X_train)

    # Implement and fit classifier
    clf = svm.SVC(kernel = "linear", C = 1) 
    clf.fit(X_train, y_train)

    pickle.dump(clf, open(current_dir + "/models/" + database_name + '-castor.joblib', 'wb'))

    # test   
    X_test, y_test      = fe.generateXYMatrice(test, k_mers, k) # OCURRENCE MATRIX
    #X_test              = fe.maxMinNormalization(X_test)
    X_test             = scaler1.transform(X_test)

    y_pred = clf.predict(X_test)        
   

    metrics = precision_recall_fscore_support(y_test, y_pred, average='weighted')
    acc = accuracy_score(y_test, y_pred)

    #print('metrics: acc, precision, recall, fscore ', acc, metrics)
    return acc, metrics[0], metrics[1], metrics[2]

    

###################################################################################################
###################################################################################################
    
# this function read the seq database proposed by me, from csv
def read_seq(path_database, database):
    path = path_database + '/' + database_name
    X = []
    X_train = []
    y_train = []
    X_test = []
    y_test = []
    #119 train, 29 test

    labels_train = np.loadtxt(path + '/train_labels.csv', dtype=(str,str), delimiter=',')
    labels_test = np.loadtxt(path + '/test_labels.csv', dtype=(str,str), delimiter=',')
    file_labels = np.vstack( (labels_train, labels_test) )
    labels = np.hstack( (labels_train[:,1], labels_test[:,1]) ) #hstack, xq tiene una dimension

    #integer_encoded, onehot_encoded, label_encoder = one_hot_encode(labels)
    
    for file_name, label in file_labels:
        file_name_without_ext = file_name.split('.')[0]
        seqs = SeqIO.parse(path + '/seq/' + file_name, "fasta") 
        
        for record in seqs:  
            X.append(str(record.seq.upper()))
            break #read one sequence by file.. all the files have one seq
    
    X = np.array(X)

    print(X.shape)
    X_train = X[0:len(labels_train)]
    X_test = X[len(labels_train):X.shape[0]]
    y_train = labels[0:len(labels_train)]
    y_test = labels[len(labels_train):X.shape[0]]

    return X_train, y_train, X_test, y_test, set(labels)


current_dir = os.path.dirname(os.path.abspath(__file__))
path_database = '/home/vicente/datasets/MLDSP/'
database_name = 'Primates'
path_database = sys.argv[1]
database_name = sys.argv[2]

# example: python3 train_kameris_castor.py '/home/vicente/datasets/MLDSP/' HIVGRPCG 0
# python3 train_kameris_castor.py '/home/vicente/datasets/MLDSP/' Primates  0
########################################################################################################

print("Reading dataset ...")
X_train, y_train, X_test, y_test, labels = read_seq(path_database, database_name)
print(X_train.shape, y_train.shape, X_test.shape, y_test.shape, labels)

k = 5

acc, presicion, recall, fscore = kameris(X_train, y_train, X_test, y_test, k, 0, database_name)
with open(current_dir + '/results_v4/results_v4.csv', "a") as myfile:
    #myfile.write("\n " + database_name + "_acc_kameris_k=" + str(k) + " " + str(acc_kameris))
    myfile.write("\n" + database_name +",kameris_k=" + str(k) + "," + str(acc) + ","+ str(presicion) + ","+ str(recall) + "," + str(fscore))

acc, presicion, recall, fscore = castor(X_train, y_train, X_test, y_test, k, 0, database_name)
with open(current_dir + '/results_v4/results_v4.csv', "a") as myfile:
    #myfile.write("\n " + database_name + "_acc_castor_k=" + str(k) + " " + str(acc_castor))
    myfile.write("\n" + database_name +",castor_k=" + str(k) + "," + str(acc) + ","+ str(presicion) + ","+ str(recall) + "," + str(fscore))
