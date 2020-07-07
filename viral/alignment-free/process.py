# replica losresultados del paper:
# DNA sequence similarity analysis using image texture analysis based on first-order statistics

from sklearn.model_selection import KFold 
from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt
from matplotlib import cm
import math
import os
import sys
import cv2
import numpy as np
import math
from scipy.stats import kurtosis, skew
from Bio import SeqIO
import pandas as pd
import seaborn as sns

from descriptor import get_features


sequences = [   'J01859.fna',       'NR_037066.fna',        'NR_040849.fna',        'NR_117152.fna',    'NR_132306.fna', 
                'NR_134817.fna',    'NR_134818.fna',        'NR_136784.fna',        'NR_148244.fna',    'NR_148787.fna', 
                'NR_152063.fna',    'KP317497.fna',         'NR_156072.fna' ]

names    = [    'Escherichia coli', 'T.Thermophilus',       'B.Wakoensis',          'T.Filiformis',     'T.Tengchongensis', 
                'S.Cameli',         'S.Tangierensis',       'T.amyloliquefaciens',  'B.Xiamenensis',    'B.Australimaris', 
                'S.Halotolerans',   'B.Maritimus',          'S.Himalayensis']

data_features = []

for sequence_file in sequences:
    data = []    
    fa_file = "sample_genomes/" + sequence_file
    sequences = SeqIO.parse(fa_file, "fasta")
    for record in sequences:
        data.append(record.seq.upper())
    seq = data[0]

    skewness, my_kurtosis, energy, entropy = get_features(seq)
    data_features.append( [skewness, my_kurtosis, energy, entropy] )
    #print([skewness, my_kurtosis, energy, entropy])


data_features = np.array(data_features)
#print(data_features)
mean_seq = data_features[0]
data_features = data_features[1:data_features.shape[0]]

#print(data_features)

distances = []
for features in data_features:
    dist = np.sqrt(np.sum((mean_seq - features)**2))
    distances.append(dist)


#sns.barplot(x=np.array(range(0,12)), y=distances)
#plt.show()

data = pd.DataFrame(data={'x': np.array(range(0,12)), 'y': distances})
grid = sns.lineplot(x=np.array(range(0,12)), y=distances)
plt.show()