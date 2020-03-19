# this script train the dataset using:
# An open-source k-mer based machine learning tool for fast and accurate [2]

import matplotlib.pyplot as plt
from sklearn.metrics import f1_score
from sklearn.feature_selection import RFE
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_recall_fscore_support
from sklearn.preprocessing import StandardScaler

from sklearn.decomposition import TruncatedSVD
from sklearn.svm import SVC
from sklearn.decomposition import PCA
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
import os
import joblib

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

###################################################################################################
#############################               kameris                 ###############################
# this funtion train kameris models proposed by
# An open-source k-mer based machine learning tool for fast and accurate [2]
def kameris(dataset, train, k, dimention_reduction):

    #minMaxScaler = MinMaxScaler(feature_range=(0, 1), copy = False)
    #scaler1 = StandardScaler()    

    #############################################
    # compute k-mer frecuences 
    X_train = []
    y_train = []
    for sample in train:
        y_train.append(sample[2])
        seq = sample[1]
        k_mers_frecuencies = kam.cgr(seq, k)   
        X_train.append(k_mers_frecuencies) 
        
    X_train = np.matrix(X_train)

    #############################################
    # scaling  
    #scaler1.fit(X_train)

    # saving scaler
    #filename = current_dir + '/models/kameris_' + dataset.split('/')[1] + "_scaler.sav"
    #joblib.dump(scaler1, filename)
    #X_train = scaler1.transform(X_train)    
    
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
    else:
        number_features = pow(4, k)

    #############################################
    # train  
    clf = SVC(kernel = "linear", C = 1) 
    clf.fit(X_train, y_train)

    return clf, number_features



#dataset_path = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/"
dataset_path = sys.argv[1]
dataset_type = sys.argv[2]
dimention_reduction = int(sys.argv[3])
k = int(sys.argv[4])

if dataset_type  == 'POL':
    datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 
                'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT'] 
elif dataset_type  == 'HIV':
    datasets = ['HIV/HIVGRPCG', 'HIV/HIVSUBCG', 'HIV/HIVSUBPOL'] 
elif dataset_type  == 'DENGE':
    datasets = ['DENGE/DENSPECG'] 

else:
    datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 
                'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT', 'HIV/HIVGRPCG', 'HIV/HIVSUBCG', 'HIV/HIVSUBPOL']

current_dir = os.path.dirname(os.path.abspath(__file__))

clf = SVC(kernel = "linear", C = 1)  

for i, dataset in enumerate(datasets):

    print(i, "\n\nEVALUATING DATASET: ", dataset, "...")
    print('===========================================================')
    print('===========================================================') 
     
    print("\n\nEvaluating with k-mer:", k, " ==========================")

    data = fe.generateLabeledData(dataset_path + dataset + "/data.fa", dataset_path  + dataset + "/class.csv")         
    
    clf, number_features = kameris(dataset, data, k, dimention_reduction)   

    filename = current_dir + '/models/kameris_' + dataset.split('/')[1] + '_dr=' + str(dimention_reduction) + "_nf=" + str(number_features) + "_k=" + str(k) + ".sav"
    print("saving model: ", filename)
    
    joblib.dump(clf, filename)


    
        
    

