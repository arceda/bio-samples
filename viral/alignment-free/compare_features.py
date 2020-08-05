# este script comprar diferente metodos de base2number

from sklearn.model_selection import KFold 
from sklearn.model_selection import train_test_split
#from matplotlib import pyplot as plt
#from matplotlib import cm
import matplotlib.pyplot as plt 
from matplotlib import pyplot
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
from descriptor import get_features_glcm
from descriptor import get_features_lbp
from descriptor import get_features_mlbp

from ete3 import PhyloTree, TreeStyle
from ete3 import Tree

from skbio import DistanceMatrix
from skbio.tree import nj


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
results_file    = current_dir + "/results/compare_features/db1"

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
results_file    = current_dir + "/results/compare_features/db2"

###################################################################################################################################
###################################################################################################################################

sequences = [   'V00662.fna',       'D38116.fna',           'D38113.fna',           'D38114.fna',       'D38115.fna', 
                'X99256.fna',       'Y18001.fna',           'X79547.fna',           'Y07726.fna',       'X63726.fna', 
                'X72004.fna',       'U20753.fna',           'X61145.fna',           'X72204.fna',       'V00654.fna', 
                'X14848.fna',       'V00711.fna',           'X83427.fna']

names    = [    'Human',            'Pygmy chimpanzee',     'Common chimpanzee',    'Gorilla',              'Orangutan', 
                'Gibbon',           'Baboon',               'Horse',                'White rhinoceros',     'Harbor seal', 
                'Gray seal',        'Cat',                  'Fin whale',            'Blue whale',           'Cow', 
                'Rat',              'Mouse',                'Platypus']

csv_mega        = current_dir + "/sample_genomes/seqs_db3_distances.csv"
seq_file_full   = current_dir + "/sample_genomes/seqs_db3.fasta"
results_file    = current_dir + "/results/compare_features/db3"


###################################################################################################################################
###################################################################################################################################

data_features_fos = []
data_features_glcm = []
data_features_lbp = []
data_features_mlbp = []

mapping_function_size = 6 # trere is 6 types of mapping functions

f_out = open(seq_file_full, "w")

for sequence_file in sequences:

    f_in = open(current_dir + "/sample_genomes/" + sequence_file, "r")
    f_out.write(f_in.read())
    f_in.close()

    data = []    
    fa_file = current_dir + "/sample_genomes/" + sequence_file
    seqs = SeqIO.parse(fa_file, "fasta")
    for record in seqs:
        data.append(record.seq.upper())      

    seq = data[0]   

    temp_fos = []
    temp_glcm = []
    temp_lbp = []
    temp_mlbp = []
    # here we evaluate each mapping funciton
    for mapping_type in range(mapping_function_size):
        skewness, my_kurtosis, energy, entropy = get_features(seq, mapping_type)
        temp_fos.append( [skewness, my_kurtosis, energy, entropy] )
        #rint("fos mapping=",mapping_type,  [skewness, my_kurtosis, energy, entropy])

        entropy, contrast, energy, correlation, homogeneity = get_features_glcm(seq, mapping_type)
        temp_glcm.append( [entropy, contrast, energy, correlation, homogeneity] )
        #print("glcm mapping=",mapping_type,  [entropy, contrast, energy, correlation, homogeneity])

        hist_lbp = get_features_lbp(seq, mapping_type)
        temp_lbp.append( hist_lbp )
        #print("lbp mapping=",mapping_type,  hist_lbp)

        hist_mlbp = get_features_mlbp(seq, mapping_type)
        temp_mlbp.append( hist_mlbp )
        #print("mlbp mapping=",mapping_type,  hist_mlbp)

    data_features_fos.append(temp_fos)
    data_features_glcm.append(temp_glcm)
    data_features_lbp.append(temp_lbp)
    data_features_mlbp.append(temp_mlbp)

f_out.close()

data_features_fos = np.array(data_features_fos)
data_features_glcm = np.array(data_features_glcm)
data_features_lbp = np.array(data_features_lbp)
data_features_mlbp = np.array(data_features_mlbp)

###################################################################################################################3
# procesamos las distancias con FOS
###################################################################################################################
full_distances_fos = []
for mapping_type in range(mapping_function_size):

    DIST_fos = np.zeros((data_features_fos.shape[0], data_features_fos.shape[0]))
    for i in range(data_features_fos.shape[0]):
        row = np.zeros(data_features_fos.shape[0])
        for j in range(i, data_features_fos.shape[0]):
            dist = np.sqrt(np.sum((data_features_fos[i][mapping_type] - data_features_fos[j][mapping_type])**2))
            row[j] = dist    
        DIST_fos[i] = row

    DIST_fos = DIST_fos + DIST_fos.T - np.diag(np.diag(DIST_fos))
    DIST_fos = (DIST_fos - np.min(DIST_fos)) / (np.max(DIST_fos) - np.min(DIST_fos))
    full_distances_fos.append( DIST_fos[0,1:DIST_fos.shape[0]] )

full_distances_fos = np.array(full_distances_fos)
print("full_distances_fos", full_distances_fos.shape)

###################################################################################################################3
# procesamos las distancias con GLCM
###################################################################################################################
full_distances_glcm = []
for mapping_type in range(mapping_function_size):

    DIST_glcm = np.zeros((data_features_glcm.shape[0], data_features_glcm.shape[0]))
    for i in range(data_features_glcm.shape[0]):
        row = np.zeros(data_features_glcm.shape[0])
        for j in range(i, data_features_glcm.shape[0]):
            dist = np.sqrt(np.sum((data_features_glcm[i][mapping_type] - data_features_glcm[j][mapping_type])**2))
            row[j] = dist    
        DIST_glcm[i] = row

    DIST_glcm = DIST_glcm + DIST_glcm.T - np.diag(np.diag(DIST_glcm))
    DIST_glcm = (DIST_glcm - np.min(DIST_glcm)) / (np.max(DIST_glcm) - np.min(DIST_glcm))
    full_distances_glcm.append( DIST_glcm[0,1:DIST_glcm.shape[0]] )

full_distances_glcm = np.array(full_distances_glcm)
print("full_distances_glcm", full_distances_glcm.shape)

###################################################################################################################3
# procesamos las distancias con LBP
###################################################################################################################
full_distances_lbp = []
for mapping_type in range(mapping_function_size):

    DIST_lbp = np.zeros((data_features_lbp.shape[0], data_features_lbp.shape[0]))
    for i in range(data_features_lbp.shape[0]):
        row = np.zeros(data_features_lbp.shape[0])
        for j in range(i, data_features_lbp.shape[0]):
            dist = np.sqrt(np.sum((data_features_lbp[i][mapping_type] - data_features_lbp[j][mapping_type])**2))
            row[j] = dist    
        DIST_lbp[i] = row

    DIST_lbp = DIST_lbp + DIST_lbp.T - np.diag(np.diag(DIST_lbp))
    DIST_lbp = (DIST_lbp - np.min(DIST_lbp)) / (np.max(DIST_lbp) - np.min(DIST_lbp))
    full_distances_lbp.append( DIST_lbp[0,1:DIST_lbp.shape[0]] )

full_distances_lbp = np.array(full_distances_lbp)
print("full_distances_lbp", full_distances_lbp.shape)

###################################################################################################################3
# procesamos las distancias con MLBP
###################################################################################################################
full_distances_mlbp = []
for mapping_type in range(mapping_function_size):

    DIST_mlbp = np.zeros((data_features_mlbp.shape[0], data_features_mlbp.shape[0]))
    for i in range(data_features_mlbp.shape[0]):
        row = np.zeros(data_features_mlbp.shape[0])
        for j in range(i, data_features_mlbp.shape[0]):
            dist = np.sqrt(np.sum((data_features_mlbp[i][mapping_type] - data_features_mlbp[j][mapping_type])**2))
            row[j] = dist    
        DIST_mlbp[i] = row

    DIST_mlbp = DIST_mlbp + DIST_mlbp.T - np.diag(np.diag(DIST_mlbp))
    DIST_mlbp = (DIST_mlbp - np.min(DIST_mlbp)) / (np.max(DIST_mlbp) - np.min(DIST_mlbp))
    full_distances_mlbp.append( DIST_mlbp[0,1:DIST_mlbp.shape[0]] )

full_distances_mlbp = np.array(full_distances_mlbp)
print("full_distances_mlbp", full_distances_mlbp.shape)

###################################################################################################################
###                    distances from mega              ###########################################################
###################################################################################################################
mega_dist_csv = pd.read_csv(csv_mega)  
mega_dist_csv = mega_dist_csv.set_index(mega_dist_csv.columns[0])
DIST_mega = mega_dist_csv.values
DIST_mega[np.isnan(DIST_mega)] = 0 # lllenamos con ceros los valores nan
DIST_mega = DIST_mega + DIST_mega.T #copiamos el triangulo inferior al superir en la matriz
distances_mega = DIST_mega[0,1:DIST_mega.shape[0]]

distances_mega = (distances_mega - np.min(distances_mega)) / (np.max(distances_mega) - np.min(distances_mega))
###################################################################################################################
###################################################################################################################

names_temp = np.array(sequences)
names_temp = names_temp[1:names_temp.shape[0]] # eliminamos el primer elemento

###################################################################################################################3
# procesamos las distancias con FOS
###################################################################################################################

plt.clf()
fig, axs = plt.subplots(3,2)
axs[0,0].plot(names_temp, full_distances_fos[0], 'b--', label='FOS-MAP0')
axs[0,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[0,0].legend(loc='upper right', fontsize=6)
axs[0,1].plot(names_temp, full_distances_fos[1], 'b--', label='FOS-MAP1')
axs[0,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[0,1].legend(loc='upper right', fontsize=6)
axs[1,0].plot(names_temp, full_distances_fos[2], 'b--', label='FOS-MAP2')
axs[1,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[1,0].legend(loc='upper right', fontsize=6)
axs[1,1].plot(names_temp, full_distances_fos[3], 'b--', label='FOS-MAP3')
axs[1,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[1,1].legend(loc='upper right', fontsize=6)
axs[2,0].plot(names_temp, full_distances_fos[4], 'b--', label='FOS-MAP4')
axs[2,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[2,0].legend(loc='upper right', fontsize=6)
axs[2,1].plot(names_temp, full_distances_fos[5], 'b--', label='FOS-MAP5')
axs[2,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[2,1].legend(loc='upper right', fontsize=6)

for ax in axs.flat:
    ax.label_outer()
    ax.yaxis.set_tick_params(labelsize=6)
    plt.sca(ax)
    plt.xticks(rotation=45, horizontalalignment='right', fontweight='light', fontsize=6 )
    plt.xlabel('Sequences', fontsize=6)

fig.text(0.04, 0.5, 'Distances', va='center', rotation='vertical', fontsize=6 )
plt.savefig( results_file + "_fos.png", dpi = 200, bbox_inches='tight')

###################################################################################################################3
# procesamos las distancias con GLCM
###################################################################################################################
plt.clf()
fig, axs = plt.subplots(3,2)
axs[0,0].plot(names_temp, full_distances_glcm[0], 'b--', label='GLCM-MAP0')
axs[0,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[0,0].legend(loc='upper right', fontsize=6)
axs[0,1].plot(names_temp, full_distances_glcm[1], 'b--', label='GLCM-MAP1')
axs[0,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[0,1].legend(loc='upper right', fontsize=6)
axs[1,0].plot(names_temp, full_distances_glcm[2], 'b--', label='GLCM-MAP2')
axs[1,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[1,0].legend(loc='upper right', fontsize=6)
axs[1,1].plot(names_temp, full_distances_glcm[3], 'b--', label='GLCM-MAP3')
axs[1,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[1,1].legend(loc='upper right', fontsize=6)
axs[2,0].plot(names_temp, full_distances_glcm[4], 'b--', label='GLCM-MAP4')
axs[2,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[2,0].legend(loc='upper right', fontsize=6)
axs[2,1].plot(names_temp, full_distances_glcm[5], 'b--', label='GLCM-MAP5')
axs[2,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[2,1].legend(loc='upper right', fontsize=6)

for ax in axs.flat:
    ax.label_outer()
    ax.yaxis.set_tick_params(labelsize=6)
    plt.sca(ax)
    plt.xticks(rotation=45, horizontalalignment='right', fontweight='light', fontsize=6 )
    plt.xlabel('Sequences', fontsize=6)

fig.text(0.04, 0.5, 'Distances', va='center', rotation='vertical', fontsize=6 )
plt.savefig( results_file + "_glcm.png", dpi = 200, bbox_inches='tight')

###################################################################################################################3
# procesamos las distancias con LBP
###################################################################################################################
plt.clf()
fig, axs = plt.subplots(3,2)
axs[0,0].plot(names_temp, full_distances_lbp[0], 'b--', label='LBP-MAP0')
axs[0,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[0,0].legend(loc='upper right', fontsize=6)
axs[0,1].plot(names_temp, full_distances_lbp[1], 'b--', label='LBP-MAP1')
axs[0,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[0,1].legend(loc='upper right', fontsize=6)
axs[1,0].plot(names_temp, full_distances_lbp[2], 'b--', label='LBP-MAP2')
axs[1,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[1,0].legend(loc='upper right', fontsize=6)
axs[1,1].plot(names_temp, full_distances_lbp[3], 'b--', label='LBP-MAP3')
axs[1,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[1,1].legend(loc='upper right', fontsize=6)
axs[2,0].plot(names_temp, full_distances_lbp[4], 'b--', label='LBP-MAP4')
axs[2,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[2,0].legend(loc='upper right', fontsize=6)
axs[2,1].plot(names_temp, full_distances_lbp[5], 'b--', label='LBP-MAP5')
axs[2,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[2,1].legend(loc='upper right', fontsize=6)

for ax in axs.flat:
    ax.label_outer()
    ax.yaxis.set_tick_params(labelsize=6)
    plt.sca(ax)
    plt.xticks(rotation=45, horizontalalignment='right', fontweight='light', fontsize=6 )
    plt.xlabel('Sequences', fontsize=6)

fig.text(0.04, 0.5, 'Distances', va='center', rotation='vertical', fontsize=6 )
plt.savefig( results_file + "_lbp.png", dpi = 200, bbox_inches='tight')

###################################################################################################################3
# procesamos las distancias con MLBP
###################################################################################################################
plt.clf()
fig, axs = plt.subplots(3,2)
axs[0,0].plot(names_temp, full_distances_mlbp[0], 'b--', label='MLBP-MAP0')
axs[0,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[0,0].legend(loc='upper right', fontsize=6)
axs[0,1].plot(names_temp, full_distances_mlbp[1], 'b--', label='MLBP-MAP1')
axs[0,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[0,1].legend(loc='upper right', fontsize=6)
axs[1,0].plot(names_temp, full_distances_mlbp[2], 'b--', label='MLBP-MAP2')
axs[1,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[1,0].legend(loc='upper right', fontsize=6)
axs[1,1].plot(names_temp, full_distances_mlbp[3], 'b--', label='MLBP-MAP3')
axs[1,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[1,1].legend(loc='upper right', fontsize=6)
axs[2,0].plot(names_temp, full_distances_mlbp[4], 'b--', label='MLBP-MAP4')
axs[2,0].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[2,0].legend(loc='upper right', fontsize=6)
axs[2,1].plot(names_temp, full_distances_mlbp[5], 'b--', label='MLBP-MAP5')
axs[2,1].plot(names_temp, distances_mega, 'r-.', label='CLUSTALW')
axs[2,1].legend(loc='upper right', fontsize=6)

for ax in axs.flat:
    ax.label_outer()
    ax.yaxis.set_tick_params(labelsize=6)
    plt.sca(ax)
    plt.xticks(rotation=45, horizontalalignment='right', fontweight='light', fontsize=6 )
    plt.xlabel('Sequences', fontsize=6)

fig.text(0.04, 0.5, 'Distances', va='center', rotation='vertical', fontsize=6 )
plt.savefig( results_file + "_mlbp.png", dpi = 200, bbox_inches='tight')


data_csv = []
error_fos = [] # save the error for each mappoing function with FOS
error_glcm = [] # save the error for each mappoing function with GLCM
error_lbp = [] # save the error for each mappoing function with LBP
error_mlbp = [] # save the error for each mappoing function with MLBP
for mapping_type in range(mapping_function_size):
    error_fos.append((np.sum((full_distances_fos[mapping_type] - distances_mega)**2))/distances_mega.shape[0])
    error_glcm.append((np.sum((full_distances_glcm[mapping_type] - distances_mega)**2))/distances_mega.shape[0])
    error_lbp.append((np.sum((full_distances_lbp[mapping_type] - distances_mega)**2))/distances_mega.shape[0])
    error_mlbp.append((np.sum((full_distances_mlbp[mapping_type] - distances_mega)**2))/distances_mega.shape[0])

data_csv.append(error_fos)
data_csv.append(error_glcm)
data_csv.append(error_lbp)
data_csv.append(error_mlbp)

data_csv = np.array(data_csv)
df = pd.DataFrame(data=data_csv.T, index=["map0", "map1", "map2", "map3", "map4", "map5"], columns=["FOS", "GLCM", "LBP", "MLBP"])
print(df)
df.to_csv(results_file + ".csv", index=True)
#print(error_fos)
#print(error_glcm)
#print(error_lbp)
#print(error_mlbp)

