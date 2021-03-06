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
from keras.layers import Conv2D, MaxPooling2D, Dense, Flatten, Dropout
from keras.utils import to_categorical
from keras.utils.vis_utils import plot_model
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import accuracy_score

from keras.utils import to_categorical
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
epochs = int(sys.argv[3])
batch_size = int(sys.argv[4])
model_type = sys.argv[5]

# example: python3 train_deep.py '/home/vicente/datasets/MLDSP/' Primates 10 32 tiny 

X_train, y_train, X_test, y_test, labels = read_cgr(5, path_database, database_name)

#create model
#epochs = 10
#batch_size = 64

##########################################################################################
##########################################################################################
if model_type == "tiny":

    model = Sequential()#add model layers
    model.add(Conv2D(32, kernel_size=3, activation='relu', input_shape=(32,32,3)))
    model.add(Conv2D(64, kernel_size=3, activation='relu'))
    model.add(Conv2D(64, kernel_size=3, activation='relu'))
    model.add(Flatten())
    model.add(Dense(len(labels), activation='softmax'))

    #compile model using accuracy to measure model performance
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

##########################################################################################
##########################################################################################


# cnn es vencido en Plants, Insects, HIV,  INFSUBMP

##########################################################################################
##########################################################################################

if model_type == "medium":

    # Initialising the CNN
    model = Sequential()

    model.add(Conv2D(input_shape=(32,32,3),filters=64,kernel_size=(3,3),padding="same", activation="relu"))

    model.add(Conv2D(filters=128,kernel_size=(3,3),padding="same", activation="relu"))

    model.add(Conv2D(filters=256,kernel_size=(3,3),padding="same", activation="relu"))

    # Step 3 - Flattening
    model.add(Flatten())
    model.add(Dense(units=64,activation="relu"))
    model.add(Dense(units = len(labels), activation = 'softmax'))

    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################

if model_type == "complex":
    # Initialising the CNN
    model = Sequential()

    model.add(Conv2D(input_shape=(32,32,3),filters=8,kernel_size=(3,3),padding="same", activation="relu"))
    #model.add(MaxPooling2D(pool_size=(2), strides=(2), padding='valid'))
    model.add(BatchNormalization())

    model.add(Conv2D(input_shape=(32,32,3),filters=16,kernel_size=(3,3),padding="same", activation="relu"))
    #model.add(MaxPooling2D(pool_size=(2), strides=(2), padding='valid'))
    model.add(BatchNormalization())

    model.add(Conv2D(input_shape=(32,32,3),filters=32,kernel_size=(3,3),padding="same", activation="relu"))
    #model.add(MaxPooling2D(pool_size=(2), strides=(2), padding='valid'))
    model.add(BatchNormalization())

    model.add(Conv2D(input_shape=(32,32,3),filters=64,kernel_size=(3,3),padding="same", activation="relu"))
    #model.add(MaxPooling2D(pool_size=(2), strides=(2), padding='valid'))
    model.add(BatchNormalization())

    model.add(Conv2D(input_shape=(32,32,3),filters=128,kernel_size=(3,3),padding="same", activation="relu"))
    #model.add(MaxPooling2D(pool_size=(2), strides=(2), padding='valid'))
    model.add(BatchNormalization())

    model.add(Flatten())

    model.add(Dense(units=256,activation="relu"))
    model.add(Dropout(0.4))
    model.add(BatchNormalization())

    model.add(Dense(units=128,activation="relu"))
    model.add(Dropout(0.4))
    model.add(BatchNormalization())

    model.add(Dense(units=64,activation="relu"))
    model.add(Dropout(0.4))
    model.add(BatchNormalization())

    model.add(Dense(units = len(labels), activation = 'softmax'))

    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

##########################################################################################
##########################################################################################
#plot_model(model, to_file='model_plot_complex.png', show_shapes=True, show_layer_names=True,  rankdir='TB', expand_nested=True, dpi=96)

#plot_model(model,  show_shapes=True, show_layer_names=True)

#sys.exit(0)

history = model.fit(X_train, y_train, epochs=epochs, batch_size=batch_size, validation_split=0.2)

model.save(current_dir + "/models/" + database_name + "-CNN-" + model_type + '.h5')


#########################################################################################################
# Plot the loss and accuracy curves for training and validation
######################################################################################################### 
plt.clf()
fig, ax = plt.subplots(2,1, figsize=(18, 10))
ax[0].plot(history.history['loss'], color='b', label="Training loss")
ax[0].plot(history.history['val_loss'], color='r', label="validation loss",axes =ax[0])
legend = ax[0].legend(loc='best', shadow=True)

ax[1].plot(history.history['accuracy'], color='b', label="Training accuracy")
ax[1].plot(history.history['val_accuracy'], color='r',label="Validation accuracy")
legend = ax[1].legend(loc='best', shadow=True)

plt.savefig(current_dir + '/results_v4/' + database_name + '_history_cnn=' + model_type + '_epoch='+ str(epochs) +'.png', dpi = 300)
#########################################################################################################
#########################################################################################################

# Confusion matrix
#plt.clf()
#fig = plt.figure(figsize=(10, 10)) # Set Figure
#y_pred = model.predict(X_test) # Predict encoded label as 2 => [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
#Y_pred = np.argmax(y_pred, 1) # Decode Predicted labels
#Y_test = np.argmax(y_test, 1) # Decode labels
# = confusion_matrix(Y_test, Y_pred) # Confusion matrix

# Plot Confusion matrix
#sns.heatmap(mat.T, square=True, annot=True, cbar=False, cmap=plt.cm.Blues)
#plt.xlabel('Predicted Values')
#plt.ylabel('True Values')
#plt.show()
#plt.savefig(current_dir + '/results/' + database_name + '_matrix_cnn=' + model_type + '_epoch='+ str(epochs) +'.png', dpi = 300)


#results = model.evaluate(X_test, y_test)
#print(results)
#print(model.metrics_names)

y_pred = model.predict(X_test)
#print(y_pred, y_test)

y_pred = np.argmax(y_pred, axis=-1)
y_test = np.argmax(y_test, axis=-1)

#print(y_pred, y_test)
acc = accuracy_score(y_test, y_pred)
metrics = precision_recall_fscore_support(y_test, y_pred, average='weighted')
#print(acc, metrics)

with open(current_dir + '/results_v4/results_v4.csv', "a") as myfile:
    myfile.write("\n" + database_name +",CNN-" + model_type + "," + str(acc) + ","+ str(metrics[0]) + ","+ str(metrics[1]) + "," + str(metrics[2]))

#with open(current_dir + '/results_v2/results_v2.txt', "a") as myfile:
#    myfile.write("\n " + database_name + "_acc_cnn=" + model_type + '_epoch='+ str(epochs) + ' ' + str(results))

#predict first 4 images in the test set
#results = model.predict(X_test)
#print(np.around(results))
#print(y_test)
