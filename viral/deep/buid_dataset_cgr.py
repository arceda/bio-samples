# this script read the dataset of MLSDP and KASTOR and generate the CRG, 
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
import csv

from sklearn.model_selection import train_test_split

# to create a SeqRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def splitMLDSP(path_database, database_name, path_dest, file_type):

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
    X = []
    y = []
        
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

            X.append(file_name)
            y.append(cluster_name)
            #if i < len(files)*0.8: # 80%                
            #    train_labels.append([file_name, cluster_name])
            #else:                
            #    test_labels.append([file_name, cluster_name])
            i += 1

    X = np.matrix(X).T   
    y = np.matrix(y).T   

      

    # con esto aseguramos que se chocoleen las muestras
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=1)

    train_labels = np.hstack((X_train, y_train))
    test_labels = np.hstack((X_test, y_test))
    np.savetxt(full_path_dest + '/train_labels.csv', train_labels, delimiter=",", fmt="%s") 
    np.savetxt(full_path_dest + '/test_labels.csv', test_labels, delimiter=",", fmt="%s") 


# esta funcion obtiene la clase de una secuencia, es el formato que usa KASTOR
def find_class(seq_id, reader_csv):
    for row in matrix_csv:
        if row[0] == seq_id:
            return row[1]

def splitKASTOR(path_database, database_name, path_dest):
    full_path_dest = path_dest + '/' + database_name
    path_seq = full_path_dest + '/' + '/seq'
    path_cgr = full_path_dest + '/' + '/cgr'

    if os.path.exists(full_path_dest) and os.path.isdir(full_path_dest):
        shutil.rmtree(full_path_dest)
    os.mkdir(full_path_dest)   

    os.mkdir(path_seq)
    os.mkdir(path_cgr)
    
    X = []
    y = []


    # read csv 
    #seqs_lables = np.loadtxt(path_database + '/' + database_name + '/class.csv', dtype=(str,str), delimiter=',')
    with open(path_database + '/' + database_name + '/class.csv') as f:
        reader = dict(csv.reader(f))


    #print(seqs_lables, seqs_lables.shape)

    # reading seqs
    sequences = SeqIO.parse(path_database + '/' + database_name + "/data.fa", "fasta")

    file_count = 0
    for record in sequences:
        seq_id  = record.id
        seq     = str(record.seq.upper()) 
        
        # create SeqRecord
        record = SeqRecord(Seq(seq, IUPAC.protein), id=seq_id, name="", description="")

        # find label
        #label = find_class(seq_id, reader_csv)
        label = reader[seq_id]

        #print("processing... ", seq_id, label)

        

        # save seq in a file
        file_name =  label + "_" + str(file_count) + ".fasta"
        output_handle = open(path_seq + "/" + file_name, "w")
        SeqIO.write(record, output_handle, "fasta")
        output_handle.close()       

        X.append(file_name)
        y.append(label)        
        
        file_count += 1

    X = np.matrix(X).T   
    y = np.matrix(y).T   
      

    # con esto aseguramos que se chocoleen las muestras
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=1)

    train_labels = np.hstack((X_train, y_train))
    test_labels = np.hstack((X_test, y_test))
    np.savetxt(full_path_dest + '/train_labels.csv', train_labels, delimiter=",", fmt="%s") 
    np.savetxt(full_path_dest + '/test_labels.csv', test_labels, delimiter=",", fmt="%s") 
   

#  this no work when ther are duplicated in seq file
'''
def splitKASTOR(path_database, database_name, path_dest):
    full_path_dest = path_dest + '/' + database_name
    path_seq = full_path_dest + '/' + '/seq'
    path_cgr = full_path_dest + '/' + '/cgr'

    if os.path.exists(full_path_dest) and os.path.isdir(full_path_dest):
        shutil.rmtree(full_path_dest)
    os.mkdir(full_path_dest)   

    os.mkdir(path_seq)
    os.mkdir(path_cgr)

    train_labels = []
    test_labels = []

    seqs_lables = np.loadtxt(path_database + '/' + database_name + '/class.csv', dtype=(str,str), delimiter=',')

    #sort by class
    seqs_lables = seqs_lables[seqs_lables[:,1].argsort()]
    types = set(seqs_lables[:,1])
    seqs_group_by_types = []
    for each_type in types:
        seqs = []
        seqs_group_by_types.append(seqs)

    #split by type or class
    for seq_id, label in seqs_lables:
        for i, each_type in enumerate(types):
            if label == each_type:
                seqs_group_by_types[i].append([seq_id, label])

    len_by_type = []
    for i, each_type in enumerate(types):
        len_by_type.append(len(seqs_group_by_types[i]))

    print(types)
    print(len_by_type)

    record_dict = SeqIO.index(path_database + '/' + database_name + "/data.fa", "fasta") # THERE IS A SEQ DUPLICATED IN influenza/INFSUBMP, and this no work

    #print(record_dict.get_raw("AJ006022").decode())
    #record_dict.close()

    file_count = 1 #use for save the sequences in different files
    for i, group in enumerate(seqs_group_by_types):
        count = 1
        for seq_id, label in group:
            #print(seq_id, label)
            #sequence = record_dict.get_raw(seq_id).decode()
            record = record_dict[seq_id]   

            # create SeqRecord
            #record = SeqRecord(Seq(sequence, IUPAC.protein), id=seq_id, name="", description="")

            #print(record)

            # save record in fasta file
            
            #file_name = str(file_count) + "_" + seq_id + ".fasta"
            file_name =  label + "_" + str(file_count) + ".fasta"
            output_handle = open(path_seq + "/" + file_name, "w")
            SeqIO.write(record, output_handle, "fasta")
            output_handle.close()

            if count < len_by_type[i]*0.8: # 80%                
                train_labels.append([file_name, label])
            else:                
                test_labels.append([file_name, label])

            count += 1
            file_count += 1


    train_labels = np.matrix(train_labels)
    test_labels = np.matrix(test_labels)
    np.savetxt(full_path_dest + '/train_labels.csv', train_labels, delimiter=",", fmt="%s") 
    np.savetxt(full_path_dest + '/test_labels.csv', test_labels, delimiter=",", fmt="%s") 


'''

def processDastaMLDSP(path_database, database_name, file_type):     
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
            #print(file_name_without_ext)      
            save_fcgr(record.id, str(record.seq.upper()), 5, path_cgr, file_name_without_ext + '_' + str(i))      
            #save_fcgr(record.id, str(record.seq.upper()), 6, path_cgr, file_name_without_ext + '_' + str(i))   
            #save_fcgr(record.id, str(record.seq.upper()), 7, path_cgr, file_name_without_ext + '_' + str(i))  
            i += 1

 
if __name__ == "__main__" :    

    current_dir = os.path.dirname(os.path.abspath(__file__))
    path_database = '/home/vicente/projects/BIOINFORMATICS/MLDSP/DataBase/'
    path_database = '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/HIV'
    path_dest = '/home/vicente/DATASETS/MLDSP/'
    database_name = 'Primates'
    path_database = sys.argv[1]
    path_dest = sys.argv[2]
    database_name = sys.argv[3]
    database_type = sys.argv[4]
    file_type = sys.argv[5]

    if database_type == 'mldsp':
        splitMLDSP(path_database, database_name, path_dest, file_type)
        processDastaMLDSP(path_dest, database_name, file_type)
    if database_type == 'kastor':
        splitKASTOR(path_database, database_name, path_dest)
        processDastaMLDSP(path_dest, database_name, file_type)

#EXAMPLES: python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' '/home/vicente/datasets/MLDSP' HIVGRPCG kastor fasta
# python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/MLDSP/DataBase/' '/home/vicente/datasets/MLDSP' Primates  mldsp  txt
