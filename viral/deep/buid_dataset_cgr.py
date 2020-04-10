# this script read the dataset of MLSDP and generate the CRG, 
# also split samples in train and test and copy thme in other folder


import collections
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import pylab
import math
import os
import shutil
import sys
import cv2
import numpy as np
from random import shuffle
from Bio import SeqIO
import re
from cgr import save_fcgr
import glob
from shutil import copyfile

def split(path_database, database_name, path_dest):

    # delete directory and crate
    full_path_dest = path_dest + '/' + database_name
    path_train = full_path_dest + '/' + '/train'
    path_test = full_path_dest + '/' + '/test'

    if os.path.exists(full_path_dest) and os.path.isdir(full_path_dest):
        shutil.rmtree(full_path_dest)
    os.mkdir(full_path_dest)   

    os.mkdir(path_train)
    os.mkdir(path_test)

    train_labels = []
    test_labels = []
        
    clusters = glob.glob(path_database + '/' + database_name + '/*' )

    for cluster in clusters: 
        cluster_name = cluster.split('/')[-1]        

        files = clusters = glob.glob(cluster + '/*.txt' )
        shuffle(files)

        i = 0
        for file in files:

            file_name = file.split('/')[-1]  
            file_name_without_ext = file_name.split('.')[0]  

            if i < len(files)*0.8: # 80%
                if os.path.exists(path_train + '/' + file_name):
                    print("archivos repetidos de diferentes clases")
                    file_name = file_name_without_ext + '_name_changed.txt' 
                copyfile(file, path_train + '/' + file_name)
                train_labels.append([file_name, cluster_name])

            else:
                if os.path.exists(path_test + '/' + file_name):
                    print("archivos repetidos de diferentes clases")
                    file_name = file_name_without_ext + '_name_changed.txt' 
                copyfile(file, path_test + '/' + file_name)
                test_labels.append([file_name, cluster_name])
            i += 1

    train_labels = np.matrix(train_labels)
    test_labels = np.matrix(train_labels)
    np.savetxt(full_path_dest + '/train_labels.csv', train_labels, delimiter=",", fmt="%s") 
    np.savetxt(full_path_dest + '/test_labels.csv', test_labels, delimiter=",", fmt="%s") 

    '''
    # delete directory and crate
    full_path_dest = path_dest + '/' + database_name
    if os.path.exists(full_path_dest) and os.path.isdir(full_path_dest):
        shutil.rmtree(full_path_dest)
    os.mkdir(full_path_dest)

    path = path_database + '/' + database_name   
    
    clusters = glob.glob(path + '/*' )

    for cluster in clusters:       

        cluster_name = cluster.split('/')[-1]

        path_train = full_path_dest + '/' + cluster_name + '/train'
        path_test = full_path_dest + '/' + cluster_name + '/test'

        os.mkdir(full_path_dest + '/' + cluster_name)
        os.mkdir(path_train)
        os.mkdir(path_test)

        files = clusters = glob.glob(cluster + '/*.txt' )

        shuffle(files)

        i = 0
        for file in files:
            file_name = file.split('/')[-1]  
            if i < len(files)*0.8: # 80%
                copyfile(file, path_train + '/' + file_name)
            else:
                copyfile(file, path_test + '/' + file_name)
            i += 1
        '''

def processDastaMLDSP(path_database, database_name): # process data in MLDSP directories
    
    path = path_database + '/' + database_name   

    number_of_clases = 0
    cluster_names = []
    points_per_cluster = []
    sequences = []
    str_all = ""
    
    #print(glob.glob(path + '/*' ))
    clusters = glob.glob(path + '/*' )
    number_of_clases = len(clusters)
    for cluster in clusters:       

        cluster_name = cluster.split('/')[-1]
        cluster_names.append( cluster_name )
       
        # read each fasta file
        files = clusters = glob.glob(cluster + '/*.txt' )
        points_per_cluster.append(len(files))
        
        # read sequences from each file, the majority have one sequence per file          
        for file in files:  
            file_name = file.split('/')[-1]  
            file_name = file_name.split('.')[0]  
            seqs = SeqIO.parse(file, "fasta") 

            i = 0          
            for record in seqs:
                
                save_fcgr(record.id, str(record.seq.upper()), 5, cluster, file_name + '_' + str(i))      
                save_fcgr(record.id, str(record.seq.upper()), 6, cluster, file_name + '_' + str(i))   
                save_fcgr(record.id, str(record.seq.upper()), 7, cluster, file_name + '_' + str(i))  
                i += 1 

    sequences_mat = np.array(sequences)

    return sequences_mat, number_of_clases, cluster_names, points_per_cluster, str_all



def processDasta(path_database, database_name): # process data in MLDSP directories
    
    path = path_database + '/' + database_name   
    
    files = glob.glob(path + '/*/*.txt' )
    for file in files:              
        file_name = file.split('/')[-1]  
        directory =  file.replace(file_name, '') 
        file_name_without_ext = file_name.split('.')[0]

        seqs = SeqIO.parse(file, "fasta") 
        i = 0          
        for record in seqs:            
            save_fcgr(record.id, str(record.seq.upper()), 5, directory, file_name_without_ext + '_' + str(i))      
            save_fcgr(record.id, str(record.seq.upper()), 6, directory, file_name_without_ext + '_' + str(i))   
            save_fcgr(record.id, str(record.seq.upper()), 7, directory, file_name_without_ext + '_' + str(i))  
            i += 1

 
if __name__ == "__main__" :    

    current_dir = os.path.dirname(os.path.abspath(__file__))
    path_database = '/home/vicente/projects/BIOINFORMATICS/MLDSP/DataBase/'
    path_dest = '/home/vicente/datasets/MLDSP/'
    database_name = 'Dengue'
    #path_database = sys.argv[1]
    #database_name = sys.argv[2]

    
    split(path_database, database_name, path_dest)
    processDasta(path_dest, database_name)