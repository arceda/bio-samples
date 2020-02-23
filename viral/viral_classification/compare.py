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

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

###################################################################################################
#############################               kameris                 ###############################
# this funtion return de acc, fcore, .... of the models proposed by
# An open-source k-mer based machine learning tool for fast and accurate [2]
def kameris(train, test, k, dimention_reduction):
    minMaxScaler = MinMaxScaler(feature_range=(0, 1), copy = False)
    scaler1 = StandardScaler()
    

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
    #X_train = minMaxScaler.fit_transform(X_train)  # this lose variance   

    scaler1.fit(X_train)
    X_train = scaler1.transform(X_train)
    #print(X_train)
    #print("STD: ", X_train.std(axis=0))
    #print("X_train > 0: ", np.count_nonzero(X_train > 0))
    #print("nfeatures: ",  np.count_nonzero(X_train > 0)*0.1  )
    #print("nfeatures: ",  (np.count_nonzero(X_train > 0)*0.1) /len(X_train) )
    #print("real nfeatues: ", X_train.shape)
    
    #############################################
    # dimention reduction   
    number_features = int((np.count_nonzero(X_train > 0)*0.1) /len(X_train) )
    if dimention_reduction == 1 and pow(4, k) > number_features and number_features > 4:
        #PCA
        #print("before PCA: ", X_train.shape)
        #pca1 = PCA(n_components = number_features, svd_solver= "arpack")
        #pca1.fit(X_train)
        #X_train = pca1.transform(X_train)
        #print("PCA X_train.shape: ", X_train.shape)

        #SVD
        svd = TruncatedSVD(n_components=number_features)
        rows, cols = X_train.shape
        svd.fit(X_train)
        X_train = svd.transform(X_train)
        print("SVD aplied ... X_train.shape: ", [rows, cols], X_train.shape)

    #############################################
    # train  
    clf = SVC(kernel = "linear", C = 1) 
    clf.fit(X_train, y_train)

    #############################################
    # compute k-mer frecuences 
    X_test = []
    y_test = []
    for sample in test:
        y_test.append(sample[2])
        seq = sample[1]
        k_mers_frecuencies = kam.cgr(seq, k)  
        X_test.append(k_mers_frecuencies)

    X_test = np.matrix(X_test)

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
    return acc, metrics[0], metrics[1], metrics[2]


###################################################################################################
#############################               castor                  ###############################
# this funtion return de acc, fcore, .... of the models proposed by
# Toward an Alignment-Free Method for Feature Extraction [1]
def castor(train, test, k, dimention_reduction, nfeatures):
    scaler1 = StandardScaler()

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
        k_mers, best_k_length = fe.getBestKmersAndFeaturesMini(train, k, range(1,15), 0.99)
        print("len_kamers_before_rfe ", len_kamers_before_rfe, " reduce fetures ", len(k_mers), k_mers)
        X_train, y_train    = fe.generateXYMatrice(train, k_mers, k)

   
        
    
    #############################################
    # scaling
    scaler1.fit(X_train)
    X_train             = scaler1.transform(X_train)

    # Implement and fit classifier
    clf = SVC(kernel = "linear", C = 1) 
    clf.fit(X_train, y_train)

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

def split_data(data, indexs):
    k_fold = []
    for i in indexs:
        k_fold.append(data[i])
    
    return k_fold

#k = 7
nfeatures = 7

#dataset_path = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/"
dataset_path = sys.argv[1]
dataset_type = sys.argv[2]
dimention_reduction = int(sys.argv[3])

if dataset_type  == 'POL':
    datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 
                'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT'] 
elif dataset_type  == 'HIV':
    datasets = ['HIV/HIVGRPCG', 'HIV/HIVSUBCG', 'HIV/HIVSUBPOL'] 
    #datasets = ['HIV/HIVSUBPOL'] 
else:
    datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 
                'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT', 'HIV/HIVGRPCG', 'HIV/HIVSUBCG', 'HIV/HIVSUBPOL']

current_dir = os.path.dirname(os.path.abspath(__file__))

clf = SVC(kernel = "linear", C = 1)  

for i, dataset in enumerate(datasets):

    print(i, "\n\nEVALUATING DATASET: ", dataset, "...")
    print('===========================================================')
    print('===========================================================') 

    csv = []

    for k in range(1,11):
        print("\n\nEvaluating with k-mer:", k, " ==========================")

        data = fe.generateLabeledData(dataset_path + dataset + "/data.fa", dataset_path  + dataset + "/class.csv")         
        
        kf = KFold(n_splits = 5, shuffle = True, random_state=1)
        i = 0

        metrics_castor = []
        metrics_kameris = []
        for train_index, test_index in kf.split(np.zeros(len(data))):            

            data_train = split_data(data, train_index)
            data_test = split_data(data, test_index)                       

            acc, pre, recall, fscore = kameris(data_train, data_test, k, dimention_reduction)
            metrics_kameris.append([acc, pre, recall, fscore])
            print("k-fold ", i, "metrics_kameris: ", acc, pre, recall, fscore)

            acc, pre, recall, fscore = castor(data_train, data_test, k, dimention_reduction, nfeatures)
            metrics_castor.append([acc, pre, recall, fscore])
            print("k-fold ", i, "metrics_castor:  ", acc, pre, recall, fscore)

            i += 1
            
        metrics_kameris         = np.matrix(metrics_kameris)
        metrics_castor          = np.matrix(metrics_castor)

        metrics_kameris_mean    = metrics_kameris.mean(0)
        metrics_castor_mean     = metrics_castor.mean(0)

        print("mean metrics_kameris: ", np.array(metrics_kameris_mean)[0])
        print("mean metrics_castor:  ", np.array(metrics_castor_mean)[0])

        row = np.append( np.array([k]), np.array(metrics_kameris_mean)[0]  )
        row = np.append( row , np.array(metrics_castor_mean)[0] )
        
        csv.append(row)

    file_name = current_dir + '/results/' + dataset.split('/')[1] + '_dr=' + str(dimention_reduction) + ".csv"
    header = "'k', 'k_acc', 'k_presicion', 'k_recall', 'k_fscore', 'c_acc', 'c_presicion', 'c_recall', 'c_fscore'"
    np.savetxt(file_name, np.array(csv), delimiter=',', fmt='%f', header=header)
    print("save file to: ", file_name)
        
    

