# this script train a CNN
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

from numpy.random import seed
seed(1)
import tensorflow as tf
tf.random.set_seed(1)



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




def numMappingPP(nucleotide):
    if nucleotide == 'A':
        return -1
    elif nucleotide == 'C':
        return 1
    elif nucleotide == 'G':
        return -1
    elif nucleotide == 'T':
        return 1


# Compute the magnitud spectrum of FFT of a numerical sequence
# seq: sequence, string of letters A, C, G, T
# median_len: median size of all sequences in a dataset
def descriptor(seq, median_len):
    #print(seq, median_len) 
    ns  = list(map(numMappingPP, seq))  # here we map the nucleotides to numbers
    ns = list(filter(None.__ne__, ns))

    ## we ensure that all the sequences have the same size 
    I   = median_len - len(ns) # change "median_len" to other length stat for length normalization
    if I > 0: # when the seq size is smaller than the 
        #print("I > 0")
        #print(ns, median_len)
        ns_temp = pywt.pad(ns, I, 'antisymmetric') #wextend('1','asym',ns,I);
        ns_temp = np.array(ns_temp)
        #print(len(ns_temp), len(ns))
        ns_new = ns_temp[I:ns_temp.shape[0]] #nsNew = nsTemp((I+1):length(nsTemp));        
    elif I < 0: # when the seq size is bigger than the median
        ns = np.array(ns)
        ns_new = ns[0:median_len]
    else:
        ns_new = np.array(ns)

    #print(ns_new.shape)
    #print("processing FFT...")
    fourier_transform = fft(ns_new)
    #print("getting  magnitude spectra...")
    magnitud_spectra = np.abs(fourier_transform) # %magnitude spectra

    return  ns_new, fourier_transform, magnitud_spectra


# read the csv and return the one hot encode of labels, use in train and test
def one_hot_encode(y):  
    # integer encode
    label_encoder = LabelEncoder()
    integer_encoded = label_encoder.fit_transform(y)
    # binary encode
    onehot_encoder = OneHotEncoder(sparse=False)
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = onehot_encoder.fit_transform(integer_encoded)

    # invert first example
    #inverted = label_encoder.inverse_transform([argmax(onehot_encoded[0, :])])
    #print(inverted)  

    return  integer_encoded, onehot_encoded, label_encoder
    
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
path_database = '/home/vicente/DATASETS/MLDSP/'
database_name = 'Primates'
path_database = sys.argv[1]
database_name = sys.argv[2]
cross_val = int(sys.argv[3])

# example: python3 train_mldsp.py '/home/vicente/DATASETS/MLDSP/' HIVGRPCG 0

print("Reading dataset ...")
X_train, y_train, X_test, y_test, labels = read_seq(path_database, database_name)
print(X_train.shape, y_train.shape, X_test.shape, y_test.shape, labels)


if cross_val == 1:
    X_train = np.vstack( (X_train, X_test) )
    y_train = np.vstack( (y_train, y_test) )


else:
    sequences_size  = list(map(len, X_train))
    median_len  = int(statistics.median(sequences_size))

    nm_val_SH     = []
    f           = []
    lg          = []

    print('Generating numerical sequences, applying DFT, computing magnitude spectra ...')

    for seq in X_train:
        ns_new, fourier_transform, magnitud_spectra = descriptor(seq, median_len)
        
        nm_val_SH.append(ns_new)
        f.append(fft(fourier_transform))
        lg.append(magnitud_spectra)
        

    #################################################################################################
    # distance calculation by Pearson correlation coefficient
    print('Computing Distance matrix .... ...')
    pearsoncorr = np.corrcoef(np.matrix(lg))
    dist_mat = (1 - pearsoncorr)/2
    print("dist_mat", dist_mat.shape)

    X_train = dist_mat

    model = svm.SVC()
    model.fit(X_train, y_train)


    ################################################TESTING##################################
    seq_distances = [] 
    for seq in X_test:
        ns_new, fourier_transform, magnitud_spectra = descriptor(seq, median_len)
        #print(magnitud_spectra, len(magnitud_spectra))
        distances = [] # the feature vector
        for observation in lg:
            #print(observation.shape)
            r = np.corrcoef(observation, magnitud_spectra)[0, 1] #corrcoef return a matrix
            distances.append((1-r)/2)

        seq_distances.append(distances)

    y_pred = model.predict(seq_distances)

    '''
    # Confusion matrix
    fig = plt.figure(figsize=(10, 10)) # Set Figure
    #print(y_pred)
    mat = confusion_matrix(y_test, y_pred) # Confusion matrix
    # Plot Confusion matrix
    sns.heatmap(mat.T, square=True, annot=True, cbar=False, cmap=plt.cm.Blues)
    plt.xlabel('Predicted Values')
    plt.ylabel('True Values')
    #plt.show()
    plt.savefig(current_dir + '/results/' + database_name + '_matrix_mldsp.png', dpi = 300)
    '''

    results = accuracy_score(y_test, y_pred)
    print(results)

    with open(current_dir + '/results/results.txt', "a") as myfile:
        myfile.write("\n " + database_name + "_acc_mldsp " + str(results))
