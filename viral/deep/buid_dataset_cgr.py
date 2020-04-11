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

def split(path_database, database_name, path_dest, file_type):

    # delete directory and crate
    full_path_dest = path_dest + '/' + database_name
    #path_train = full_path_dest + '/' + '/train'
    #path_test = full_path_dest + '/' + '/test'
    path_seq = full_path_dest + '/' + '/seq'
    path_cgr = full_path_dest + '/' + '/cgr'

    if os.path.exists(full_path_dest) and os.path.isdir(full_path_dest):
        shutil.rmtree(full_path_dest)
    os.mkdir(full_path_dest)   

    os.mkdir(path_seq)
    os.mkdir(path_cgr)

    train_labels = []
    test_labels = []
        
    clusters = glob.glob(path_database + '/' + database_name + '/*' )

    for cluster in clusters: 
        cluster_name = cluster.split('/')[-1]        

        files = clusters = glob.glob(cluster + '/*.' + file_type )
        shuffle(files)

        i = 0
        for file in files:
            file_name = file.split('/')[-1]  
            file_name_without_ext = file_name.split('.')[0]  

            if os.path.exists(path_seq + '/' + file_name):
                print("archivos repetidos de diferentes clases")
                file_name = file_name_without_ext + '_name_changed.txt' 
            copyfile(file, path_seq + '/' + file_name)

            if i < len(files)*0.8: # 80%                
                train_labels.append([file_name, cluster_name])
            else:                
                test_labels.append([file_name, cluster_name])
            i += 1

    train_labels = np.matrix(train_labels)
    test_labels = np.matrix(test_labels)
    np.savetxt(full_path_dest + '/train_labels.csv', train_labels, delimiter=",", fmt="%s") 
    np.savetxt(full_path_dest + '/test_labels.csv', test_labels, delimiter=",", fmt="%s") 


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



def processDasta(path_database, database_name, file_type):     
    path = path_database + '/' + database_name      
    path_cgr = path + '/' + '/cgr'

    if os.path.exists(path_cgr) and os.path.isdir(path_cgr):
        shutil.rmtree(path_cgr)
    os.mkdir(path_cgr)  
    
    files = glob.glob(path + '/*/*.'  + file_type )
    for file in files:              
        file_name = file.split('/')[-1]  
        directory =  file.replace(file_name, '') 
        file_name_without_ext = file_name.split('.')[0]

        seqs = SeqIO.parse(file, "fasta") 
        i = 0          
        for record in seqs:            
            save_fcgr(record.id, str(record.seq.upper()), 5, path_cgr, file_name_without_ext + '_' + str(i))      
            #save_fcgr(record.id, str(record.seq.upper()), 6, path_cgr, file_name_without_ext + '_' + str(i))   
            #save_fcgr(record.id, str(record.seq.upper()), 7, path_cgr, file_name_without_ext + '_' + str(i))  
            i += 1

 
if __name__ == "__main__" :    

    current_dir = os.path.dirname(os.path.abspath(__file__))
    path_database = '/home/vicente/projects/BIOINFORMATICS/MLDSP/DataBase/'
    path_dest = '/home/vicente/datasets/MLDSP/'
    database_name = 'Primates'
    path_database = sys.argv[1]
    path_dest = sys.argv[2]
    database_name = sys.argv[3]
    file_type = sys.argv[4]

    
    split(path_database, database_name, path_dest, file_type)
    processDasta(path_dest, database_name, file_type)