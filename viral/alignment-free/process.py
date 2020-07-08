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

current_dir = os.path.dirname(os.path.abspath(__file__))

###################################################################################################################################
###################################################################################################################################

sequences = [   'J01859.fna',       'NR_037066.fna',        'NR_040849.fna',        'NR_117152.fna',    'NR_132306.fna', 
                'NR_134817.fna',    'NR_134818.fna',        'NR_136784.fna',        'NR_148244.fna',    'NR_148787.fna', 
                'NR_152063.fna',    'KP317497.fna',         'NR_156072.fna' ]

names    = [    'Escherichia coli', 'T.Thermophilus',       'B.Wakoensis',          'T.Filiformis',     'T.Tengchongensis', 
                'S.Cameli',         'S.Tangierensis',       'T.amyloliquefaciens',  'B.Xiamenensis',    'B.Australimaris', 
                'S.Halotolerans',   'B.Maritimus',          'S.Himalayensis']

csv_mega        = current_dir + "/sample_genomes/seqs_db1_distances.csv"
seq_file_full   = current_dir + "/sample_genomes/seqs_db1.fasta"
results_file    = current_dir + "/results/db1.png"

###################################################################################################################################
###################################################################################################################################

sequences = [   'L00016.fna',       'M22650.fna',           'M22651.fna',           'M22653.fna',       'M22654.fna', 
                'M22655.fna',       'M22656.fna',           'M22657.fna',           'V00658.fna',       'V00659.fna', 
                'V00672.fna',       'V00675.fna']

names    = [    'Human',            'Macaca mulatta',       'Macaca fuscata',       'Macaca fascicularis',  'Macaca sylvanus', 
                'Saimiri sciureus', 'Tarsius syrichta',     'Lemur catta',          'Gorilla',              'Hylobates', 
                'Chimpanzee',       'Sumatran Orangutan']

csv_mega        = current_dir + "/sample_genomes/seqs_db2_distances.csv"
seq_file_full   = current_dir + "/sample_genomes/seqs_db2.fasta"
results_file    = current_dir + "/results/db2.png"

###################################################################################################################################
###################################################################################################################################

data_features = []

f_out = open(seq_file_full, "w")

for sequence_file in sequences:

    f_in = open(current_dir + "/sample_genomes/" + sequence_file, "r")
    f_out.write(f_in.read())
    f_in.close()

    data = []    
    fa_file = current_dir + "/sample_genomes/" + sequence_file
    sequences = SeqIO.parse(fa_file, "fasta")
    for record in sequences:
        #print(record.id)
        #print(record.description)
        data.append(record.seq.upper())      

    seq = data[0]    

    skewness, my_kurtosis, energy, entropy = get_features(seq)
    data_features.append( [skewness, my_kurtosis, energy, entropy] )
    #print([skewness, my_kurtosis, energy, entropy])
f_out.close()

data_features = np.array(data_features)
#print(data_features)
mean_seq = data_features[0]
data_features = data_features[1:data_features.shape[0]]

#print(data_features)

#distances = []
#for features in data_features:
#    dist = np.sqrt(np.sum((mean_seq - features)**2))
#    distances.append(dist)

distances = (distances - np.min(distances)) / (np.max(distances) - np.min(distances))

### distances from mega ###########################################################
mega_dist_csv = pd.read_csv(csv_mega)  
mega_dist_csv = mega_dist_csv.set_index(mega_dist_csv.columns[0])
dist_data = mega_dist_csv.values
dist_data = dist_data[ 1:dist_data.shape[0], 0:dist_data.shape[0]-1]
distances_mega = dist_data[:, 0]

distances_mega = (distances_mega - np.min(distances_mega)) / (np.max(distances_mega) - np.min(distances_mega))

names = np.array(names)
names = names[1:names.shape[0]] # eliminamos el primer elemento
#print(names)

plt.clf()
plt.plot(names, distances, label='Proposed method')
plt.plot(names, distances_mega, label='Mega')
#chart = sns.lineplot(x = names, y = distances, label='Proposed method')
#chart = sns.lineplot(x = names, y = distances_mega, label='Mega')
plt.xticks(rotation=45, horizontalalignment='right', fontweight='light', fontsize=6 )
plt.legend(loc='upper right')
#plt.xlabel('longitude')
#plt.ylabel('latitude')
#plt.show()
plt.savefig( results_file, dpi = 200, bbox_inches='tight')
