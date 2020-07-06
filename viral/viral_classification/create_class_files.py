#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:05:54 2020

@author: vicente
"""

# this script labels the incomplete datasets in csv files

import numpy as np
from Bio import SeqIO
import re
import glob
from iteration_utilities import duplicates


def checkIfDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        print(list(duplicates(listOfElems)))
        return True


# labels datasets
labels_dat = ['PAPILLOMA/HPVGENCG', 'PAPILLOMA/HPVSPECG', 'HIV/HIVGRPCG', 'HIV/HIVSUBCG', 'HIV/HIVSUBPOL', 'HEPATITIS-B/HBVGENCG'] # from CASTOR dataset
unlabels_dat = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT', 'RHINOVIRUS/RHISPECG', 'DENGE/DENSPECG', 'INFLUENZA/INSUBFNA', 'INFLUENZA/INFSUBHA', 'INFLUENZA/INFSUBMP', 'EBOLA/EBOSPECG']

path_dataset ='/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/'

#print(labels_dat)
#print(unlabels_dat)

for i, folder in enumerate(unlabels_dat):
    path = path_dataset + folder
    
    # testing if there are all the files
    #print(path)
    #files = glob.glob(path + '/class.csv')
    #files = glob.glob(path + '/class.csv')
    #print(files)   

        
    # PARSING POLYOMAVIRUS y Rhinovirus    
    if i == 0 or i == 1 or i == 2 or i == 3 or i == 4 or i == 5:
        # utilized for extract substring(label) from de id sequence
        flag_ini = "|" 
        flag_fin = "_"
    
    # PARSING DENGE
    elif i == 6:
        # utilized for extract substring(label) from de id sequence
        flag_ini = "Dengue_" 
        flag_fin = "|"
    
    # PARSING INFLUENZA
    elif i == 7 or i == 8 or i == 9:
        # utilized for extract substring(label) from de id sequence
        flag_ini = "|" 
        flag_fin = "|"
        
    # PARSING EBOLA
    elif i == 10:
        # utilized for extract substring(label) from de id sequence
        flag_ini = "" 
        flag_fin = "_ebola"
    
    else:
        break
        
    data = []
    labels = []
    ides = []
    
    seq_file = path_dataset + folder + '/data.fa'    
    
    print('\nparsing ... ', seq_file)
    
    sequences = SeqIO.parse(seq_file, "fasta")
    for record in sequences:
        data.append([record.id, record.seq.upper()])
        record_id = record.id
        
        #print(record_id)
        # extract label from "1|BK_VP1" , the label is BK
        # the label is between "|" and "_"        
        
        if i == 6: # DENGE
            lbl = record_id.partition(flag_ini)[2].partition(flag_fin)[0]
            lbl = lbl[0:7]            
        if i == 10: # EBOLA
            lbl = record_id.partition(flag_fin)[0]
        else:
            lbl = record_id.partition(flag_ini)[2].partition(flag_fin)[0]
            
        labels.append(lbl)        
        ides.append(record_id)
    
    # check data 
    labels_unique = set(labels)
    print("duplicate ides: ", checkIfDuplicates(ides))
    print('unique_labels', len(labels_unique), labels_unique)
    print('instances', len(ides))
    
    csv = np.column_stack((ides, labels))
    np.savetxt(path_dataset + folder + '/class.csv', csv, delimiter=',', fmt='%s') 


'''
flag_ini = "|" 
flag_fin = "|"
data = []
labels = []
ides = []

seq_file = path_dataset + unlabels_dat[9] + '/data.fa'    

print('\nparsing ... ', seq_file)

sequences = SeqIO.parse(seq_file, "fasta")
for record in sequences:
    data.append([record.id, record.seq.upper()])
    record_id = record.id
    
    #print(record_id)
    # extract label from "1|BK_VP1" , the label is BK
    # the label is between "|" and "_"
    
    lbl = record_id.partition(flag_ini)[2].partition(flag_fin)[0]
    #lbl = record_id.partition((flag_fin)[0]
    
    labels.append(lbl)
    ides.append(record_id)
    
# check data 
labels_unique = set(labels)
print("duplicate ides: ", checkIfDuplicates(ides))
print('unique_labels', len(labels_unique), labels_unique)
print('instances', len(ides))


csv = np.column_stack((ides, labels))
#np.savetxt(path_dataset + folder + '/class.csv', csv, delimiter=',', fmt='%s') 

'''