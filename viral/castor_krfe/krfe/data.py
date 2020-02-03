# Imports
import csv
import numpy as np
from Bio import SeqIO

def generateLabeledData(fa_file, cls_file):
    # Variable data 
    data = []

    # Open the class file
    with open(cls_file) as f:
        reader = dict(csv.reader(f))

    # Open the sequences file
    for record in SeqIO.parse(fa_file, "fasta"):
        if record.id in reader:
            # Generate table [Id, Sequences, Class]
            data.append([record.id, record.seq.upper(), reader[record.id]])

    # Return data
    return data


def generateData(fa_file):
    # Variable data 
    data = []

    # Open the sequences file
    for record in SeqIO.parse(fa_file, "fasta"):
        # Generate table [Id, Sequences, Class]
        data.append([record.id, record.seq.upper(), None])

    # Return data
    return data

if __name__ == "__main__" :
    #training_data = generateLabeledData("../Data/HIVGRPCG/data.fa", "../Data/HIVGRPCG/class.csv")
    
    fa_file = "../Data/HIVGRPCG/data.fa"
    cls_file = "../Data/HIVGRPCG/class.csv"
    
    data = []

    # Open the class file
    with open(cls_file) as f:
        reader = dict(csv.reader(f))

    # Open the sequences file}
    sequences = SeqIO.parse(fa_file, "fasta")
    print(len(sequences))
    for record in sequences:
        if record.id in reader:
            # Generate table [Id, Sequences, Class]
            data.append([record.id, record.seq.upper(), reader[record.id]])