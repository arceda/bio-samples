# this script train the dataset using:
# An open-source k-mer based machine learning tool for fast and accurate [2]

import sys
import numpy as np

import mykameris as kam
import os
import joblib
from Bio import SeqIO
import re

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

#fasta = "/home/vicente/projects/BIOINFORMATICS/bio-samples/viral/viral_classification/hiv1-genomes/A1.fasta"
#dataset = HIVSUBCG
fasta = sys.argv[1]
dataset = sys.argv[2]

current_dir = os.path.dirname(os.path.abspath(__file__))

#folder = "/home/vicente/projects/BIOINFORMATICS/bio-samples/viral/viral_classification/models/"
if dataset[0:3] == "HIV":
    k = 5
    model = current_dir + "/models/kameris_" + dataset + "_dr=0_nf=1024_k=5.sav"
elif dataset[0:3] == "POL":
    k = 2
    model = current_dir + "/models/kameris_" + dataset + "_dr=0_nf=16_k=2.sav"




clf = joblib.load(model)

data = []
sequences = SeqIO.parse(fasta, "fasta")
for record in sequences:
    data.append([record.id, record.seq.upper()])

X_test = []
for seq in data:
    #print("processing seq: ", seq[0])
    k_mers_frecuencies = kam.cgr(seq[1], k)  
    X_test.append(k_mers_frecuencies)

result = clf.predict(X_test)
print(result)

    
        
    

