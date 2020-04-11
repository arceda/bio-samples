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
from numpy import argmax
import matplotlib.pyplot as plt

import statistics
import pywt
import csv


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

    integer_encoded, onehot_encoded, label_encoder = one_hot_encode(labels)
    
    for file_name, label in file_labels:
        file_name_without_ext = file_name.split('.')[0]
        seqs = SeqIO.parse(path + '/seq/' + file, "fasta") 
        
        for record in seqs:  
            X.append(record.seq.upper())
            break #read one sequence by file.. all the files have one seq
    
    X = np.array(X)
    X_train = X[0:len(labels_train),:,:,:]
    X_test = X[len(labels_train):X.shape[0],:,:,:]
    y_train = onehot_encoded[0:len(labels_train),:]
    y_test = onehot_encoded[len(labels_train):X.shape[0],:]

    return X_train, y_train, X_test, y_test, set(labels)

def read_cgr(k, path_database, database_name):
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

    integer_encoded, onehot_encoded, label_encoder = one_hot_encode(labels)
    
    for file_name, label in file_labels:
        file_name_without_ext = file_name.split('.')[0]
        # por ahora solo estamos considerando una secuenica por fasta
        img = cv2.imread(path + '/cgr/' + file_name_without_ext + '_0_k=' + str(k) + '.jpg' ) 
        #cv2.imshow("win", img)
        #cv2.waitKey()
        X.append(img)

    X = np.array(X)
    X_train = X[0:len(labels_train),:,:,:]
    X_test = X[len(labels_train):X.shape[0],:,:,:]
    y_train = onehot_encoded[0:len(labels_train),:]
    y_test = onehot_encoded[len(labels_train):X.shape[0],:]
    
    #print(X.shape, X_train.shape, X_test.shape, y_train.shape, y_test.shape)
    return X_train, y_train, X_test, y_test, set(labels)

current_dir = os.path.dirname(os.path.abspath(__file__))
path_database = '/home/vicente/datasets/MLDSP/'
database_name = 'Primates'
path_database = sys.argv[1]
database_name = sys.argv[2]

X_train, y_train, X_test, y_test, labels = read_cgr(5, path_database, database_name)

#create model
epochs = 10
batch_size = 64


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
#model.fit(X_train, y_train, epochs=3)


'''
model=Sequential()
#model.add(Lambda(standardize,input_shape=(28,28,1)))    
model.add(Conv2D(filters=64, kernel_size = (3,3), activation="relu", input_shape=(32,32,3)))
model.add(Conv2D(filters=64, kernel_size = (3,3), activation="relu"))

model.add(MaxPooling2D(pool_size=(2,2)))
model.add(BatchNormalization())
model.add(Conv2D(filters=128, kernel_size = (3,3), activation="relu"))
model.add(Conv2D(filters=128, kernel_size = (3,3), activation="relu"))

model.add(MaxPooling2D(pool_size=(2,2)))
model.add(BatchNormalization())    
model.add(Conv2D(filters=256, kernel_size = (3,3), activation="relu"))
    
model.add(MaxPooling2D(pool_size=(2,2)))
    
model.add(Flatten())
model.add(BatchNormalization())
model.add(Dense(512,activation="relu"))
    
model.add(Dense(len(labels),activation="softmax"))    
model.compile(loss="categorical_crossentropy", optimizer="adam", metrics=["accuracy"])
'''

history = model.fit(X_train, y_train, epochs=epochs, validation_split=0.2)

#plot_model(model, to_file='model_plot.png', show_shapes=True, show_layer_names=True)

# Plot the loss and accuracy curves for training and validation
# ##################################################################### 
#fig, ax = plt.subplots(2,1, figsize=(18, 10))
#ax[0].plot(history.history['loss'], color='b', label="Training loss")
#ax[0].plot(history.history['val_loss'], color='r', label="validation loss",axes =ax[0])
#legend = ax[0].legend(loc='best', shadow=True)

#ax[1].plot(history.history['accuracy'], color='b', label="Training accuracy")
#ax[1].plot(history.history['val_accuracy'], color='r',label="Validation accuracy")
#legend = ax[1].legend(loc='best', shadow=True)



# Confusion matrix
fig = plt.figure(figsize=(10, 10)) # Set Figure
y_pred = model.predict(X_test) # Predict encoded label as 2 => [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
Y_pred = np.argmax(y_pred, 1) # Decode Predicted labels
Y_test = np.argmax(y_test, 1) # Decode labels
mat = confusion_matrix(Y_test, Y_pred) # Confusion matrix

# Plot Confusion matrix
sns.heatmap(mat.T, square=True, annot=True, cbar=False, cmap=plt.cm.Blues)
plt.xlabel('Predicted Values')
plt.ylabel('True Values')
#plt.show()
plt.savefig(current_dir + '/results/' + database_name + '_matrix_cnn=tiny_epoch=10.png', dpi = 300)


results = model.evaluate(X_test, y_test)
print(results)

with open(current_dir + '/results/results.txt', "a") as myfile:
    myfile.write("\n " + database_name + "_acc_cnn=tiny_epoch=10 " + str(results))

#predict first 4 images in the test set
#results = model.predict(X_test)
#print(np.around(results))
#print(y_test)