# this script train a CNN

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
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from numpy import argmax
import matplotlib.pyplot as plt

import csv

def one_hot_encode(y_train, y_test)
    # integer encode
    label_encoder = LabelEncoder()
    integer_encoded = label_encoder.fit_transform(y_train)
    print(integer_encoded)
    # binary encode
    onehot_encoder = OneHotEncoder(sparse=False)
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = onehot_encoder.fit_transform(integer_encoded)
    print(onehot_encoded)
    # invert first example
    inverted = label_encoder.inverse_transform([argmax(onehot_encoded[0, :])])
    print(inverted)   
    

def read_cgr(k, path_database, database_name):
    path = path_database + '/' + database_name
    X_train = []
    y_train = []
    X_test = []
    y_test = []

    with open(path + '/train_labels.csv') as f:
        labels_train_reader =  (csv.reader(f))
    with open(path + '/test_labels.csv') as f:
        labels_test_reader = dict(csv.reader(f))

    
    files = glob.glob(path + '/train/*.txt' )
    for file in files: 
        file_name = file.split('/')[-1]  
        directory =  file.replace(file_name, '') 
        file_name_without_ext = file_name.split('.')[0]
        print(file_name)

        img = cv2.imread(file_name_without_ext + '_0_k=' + str(k) + '.jpg' ) # por ahora solo estamos considerando una secuenica por fasta
        X_train.append(img)
        y_train.append(labels_train_reader[file_name])


    files = glob.glob(path + '/test/*.txt' )
    for file in files: 
        file_name = file.split('/')[-1]  
        directory =  file.replace(file_name, '') 
        file_name_without_ext = file_name.split('.')[0]
        print(file_name)

        img = cv2.imread(file_name_without_ext + '_0_k=' + str(k) + '.jpg' ) # por ahora solo estamos considerando una secuenica por fasta
        X_test.append(img)
        y_test.append(labels_test_reader[file_name])



    return np.array(X_train), np.array(y_train), np.array(X_test), np.array(y_test)

current_dir = os.path.dirname(os.path.abspath(__file__))
path_database = '/home/vicente/datasets/MLDSP/'
database_name = 'Primates'
#path_database = sys.argv[1]
#database_name = sys.argv[2]

X_train, y_train = read_cgr(5, path_database, database_name)
print(X_train.shape, y_train.shape)
print(y_train)



#create model
model = Sequential()#add model layers
model.add(Conv2D(32, kernel_size=3, activation='relu', input_shape=(32,32,3)))
model.add(Conv2D(64, kernel_size=3, activation='relu'))
model.add(Conv2D(64, kernel_size=3, activation='relu'))
model.add(Flatten())
model.add(Dense(len(labels), activation='softmax'))


#compile model using accuracy to measure model performance
model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

#train the model
#model.fit(X_train, y_train, validation_data=(X_test, y_test), epochs=3)
model.fit(X_train, y_train, epochs=3)

#predict first 4 images in the test set
results = model.predict(X_test[:4])
print(round(results))
print(y_test)