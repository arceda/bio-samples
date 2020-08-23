import os
import subprocess
import numpy as np
from Bio import SeqIO
import re

current_dir = os.path.dirname(os.path.abspath(__file__))

def read_fasta(fa_file): 
    data = []

    sequences = SeqIO.parse(fa_file, "fasta")

    for record in sequences:
        data.append([record.id, record.seq.upper()])
        #print(record.id)
        #print(record.seq.upper())     

    return data

def cgr(seq, k, order='ACGT'):
    result = np.zeros(pow(4, k), dtype='uint8')
    #print(result.shape, result)

    x = pow(2, k-1)
    y = pow(2, k-1)

    #print("len seq:", len(seq))
    for i, c  in enumerate(seq):
    #for i, c  in enumerate(range(20)):
        #print(i, c)
        if seq[i] != order[0] and seq[i] != order[1] and seq[i] != order[2] and seq[i] != order[3]:
            continue

        x >>= 1  #x /= 2
        if seq[i] == order[2] or seq[i] == order[3]:
            x += pow(2, k-1)

        y >>= 1  #y /= 2
        if seq[i] == order[0] or seq[i] == order[3]:
            y += pow(2, k-1)
        
        
        if i >= (k - 1):
            #print((pow(2, k))*y + x)
            pos = int( pow(2, k)*y + x )
            #print(pos)
            result[pos] += 1

    return result


