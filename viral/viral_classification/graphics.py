import numpy as np
import sys
import csv
import matplotlib as mpl 
## agg backend is used to create plot as a .png file
#mpl.use('agg')
import matplotlib.pyplot as plt 
from matplotlib import pyplot
import os
import pandas as pd

current_dir = os.path.dirname(os.path.abspath(__file__))
#datasets = ['POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST', 'POLSPELT'] 
datasets = ['HIVGRPCG', 'HIVSUBCG', 'HIVSUBPOL'] 

datasets = ['POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST', 'POLSPELT', 'HIVGRPCG', 'HIVSUBCG', 'HIVSUBPOL'] 


#dataset_path = sys.argv[1]
#dataset_path = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/"

# each graphic separated
############################################################################################
'''
for i, dataset in enumerate(datasets):
    print("\n\n", i, "EVALUATING DATASET: ", dataset, "...\n")    

    # READ WITH PANDAS
    #df = pd.read_csv(file_name, skiprows=1)
    #print(df.head())
    #print(df)

    file_name = current_dir + '/results/' + dataset + '_dr=0.csv'
    data = np.genfromtxt(file_name, delimiter=',')

    #print(data)

    plt.plot(data[:, 0], data[:, 4], 'r', label='kameris')
    plt.plot(data[:, 0], data[:, 8], 'b--', label='castor')
    #plt.axis('equal')
    plt.xlabel('k')
    plt.ylabel('fscore')
    plt.legend()
    plt.show()
'''
############################################################################################


POLSPEVP1 = np.genfromtxt(current_dir + '/results/POLSPEVP1_dr=0.csv', delimiter=',')
POLSPEVP2 = np.genfromtxt(current_dir + '/results/POLSPEVP2_dr=0.csv', delimiter=',')
POLSPEVP3 = np.genfromtxt(current_dir + '/results/POLSPEVP3_dr=0.csv', delimiter=',')
POLSPEST = np.genfromtxt(current_dir + '/results/POLSPEST_dr=0.csv', delimiter=',')

POLSPELT = np.genfromtxt(current_dir + '/results/POLSPELT_dr=0.csv', delimiter=',')
HIVGRPCG = np.genfromtxt(current_dir + '/results/HIVGRPCG_dr=0.csv', delimiter=',')
HIVSUBCG = np.genfromtxt(current_dir + '/results/HIVSUBCG_dr=0.csv', delimiter=',')
HIVSUBPOL = np.genfromtxt(current_dir + '/results/HIVSUBPOL_dr=0.csv', delimiter=',')

POLSPEVP1_dr = np.genfromtxt(current_dir + '/results/POLSPEVP1_dr=1.csv', delimiter=',')
POLSPEVP2_dr = np.genfromtxt(current_dir + '/results/POLSPEVP2_dr=1.csv', delimiter=',')
POLSPEVP3_dr = np.genfromtxt(current_dir + '/results/POLSPEVP3_dr=1.csv', delimiter=',')
POLSPEST_dr = np.genfromtxt(current_dir + '/results/POLSPEST_dr=1.csv', delimiter=',')

POLSPELT_dr = np.genfromtxt(current_dir + '/results/POLSPELT_dr=1.csv', delimiter=',')

############################################################################################
# comparison between 'POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST' without dimentionality reduction
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(POLSPEVP1[:, 0], POLSPEVP1[:, 4], 'r', label='kameris')
axs[0, 0].plot(POLSPEVP1[:, 0], POLSPEVP1[:, 8], 'b--', label='castor')
axs[0, 0].legend(loc='lower right')
axs[0, 0].set_title('POLSPEVP1')
axs[0, 1].plot(POLSPEVP2[:, 0], POLSPEVP2[:, 4], 'r', label='kameris')
axs[0, 1].plot(POLSPEVP2[:, 0], POLSPEVP2[:, 8], 'b--', label='castor')
axs[0, 1].legend(loc='lower right')
axs[0, 1].set_title('POLSPEVP2')
axs[1, 0].plot(POLSPEVP3[:, 0], POLSPEVP3[:, 4], 'r', label='kameris')
axs[1, 0].plot(POLSPEVP3[:, 0], POLSPEVP3[:, 8], 'b--', label='castor')
axs[1, 0].legend(loc='lower right')
axs[1, 0].set_title('POLSPEVP3')
axs[1, 1].plot(POLSPEST[:, 0], POLSPEST[:, 4], 'r', label='kameris')
axs[1, 1].plot(POLSPEST[:, 0], POLSPEST[:, 8], 'b--', label='castor')
axs[1, 1].legend(loc='lower right')
axs[1, 1].set_title('POLSPEST')

for ax in axs.flat:
    ax.set(xlabel='k', ylabel='fscore')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

#plt.show()
plt.savefig(current_dir + '/results/' + 'comparison1_kameris_castor_dr=0.png')
#############################################################################################


############################################################################################
# comparison between 'POLSPELT', 'HIVGRPCG', 'HIVSUBCG', 'HIVSUBPOL' without dimentionality reduction
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(POLSPELT[:, 0], POLSPELT[:, 4], 'r', label='kameris')
axs[0, 0].plot(POLSPELT[:, 0], POLSPELT[:, 8], 'b--', label='castor')
axs[0, 0].legend(loc='lower right')
axs[0, 0].set_title('POLSPELT')
axs[0, 1].plot(HIVGRPCG[:, 0], HIVGRPCG[:, 4], 'r', label='kameris')
axs[0, 1].plot(HIVGRPCG[:, 0], HIVGRPCG[:, 8], 'b--', label='castor')
axs[0, 1].legend(loc='lower right')
axs[0, 1].set_title('HIVGRPCG')
axs[1, 0].plot(HIVSUBCG[:, 0], HIVSUBCG[:, 4], 'r', label='kameris')
axs[1, 0].plot(HIVSUBCG[:, 0], HIVSUBCG[:, 8], 'b--', label='castor')
axs[1, 0].legend(loc='lower right')
axs[1, 0].set_title('HIVSUBCG')
axs[1, 1].plot(HIVSUBPOL[:, 0], HIVSUBPOL[:, 4], 'r', label='kameris')
axs[1, 1].plot(HIVSUBPOL[:, 0], HIVSUBPOL[:, 8], 'b--', label='castor')
axs[1, 1].legend(loc='lower right')
axs[1, 1].set_title('HIVSUBPOL')

for ax in axs.flat:
    ax.set(xlabel='k', ylabel='fscore')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

#plt.show()
plt.savefig(current_dir + '/results/' + 'comparison2_kameris_castor_dr=0.png')
#############################################################################################


############################################################################################
# comparison between kameris with and with SVD 'POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST' 
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(POLSPEVP1[:, 0], POLSPEVP1[:, 4], 'r', label='kameris')
axs[0, 0].plot(POLSPEVP1_dr[:, 0], POLSPEVP1_dr[:, 4], 'b--', label='kameris-svd')
axs[0, 0].legend(loc='lower right')
axs[0, 0].set_title('POLSPEVP1')
axs[0, 1].plot(POLSPEVP2[:, 0], POLSPEVP2[:, 4], 'r', label='kameris')
axs[0, 1].plot(POLSPEVP2_dr[:, 0], POLSPEVP2_dr[:, 4], 'b--', label='kameris-svd')
axs[0, 1].legend(loc='lower right')
axs[0, 1].set_title('POLSPEVP2')
axs[1, 0].plot(POLSPEVP3[:, 0], POLSPEVP3[:, 4], 'r', label='kameris')
axs[1, 0].plot(POLSPEVP3_dr[:, 0], POLSPEVP3_dr[:, 4], 'b--', label='kameris-svd')
axs[1, 0].legend(loc='lower right')
axs[1, 0].set_title('POLSPEVP3')
axs[1, 1].plot(POLSPEST[:, 0], POLSPEST[:, 4], 'r', label='kameris')
axs[1, 1].plot(POLSPEST_dr[:, 0], POLSPEST_dr[:, 4], 'b--', label='kameris-svd')
axs[1, 1].legend(loc='lower right')
axs[1, 1].set_title('POLSPEST')

for ax in axs.flat:
    ax.set(xlabel='k', ylabel='fscore')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

#plt.show()
plt.savefig(current_dir + '/results/' + 'comparison3_kameris_dr=0_1.png')
#############################################################################################



############################################################################################
# comparison between castor with and with SVD 'POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST' 
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(POLSPEVP1[:, 0], POLSPEVP1[:, 8], 'r', label='castor')
axs[0, 0].plot(POLSPEVP1_dr[:, 0], POLSPEVP1_dr[:, 8], 'b--', label='castor-rfe')
axs[0, 0].legend(loc='lower right')
axs[0, 0].set_title('POLSPEVP1')
axs[0, 1].plot(POLSPEVP2[:, 0], POLSPEVP2[:, 8], 'r', label='castor')
axs[0, 1].plot(POLSPEVP2_dr[:, 0], POLSPEVP2_dr[:, 8], 'b--', label='castor-rfe')
axs[0, 1].legend(loc='lower right')
axs[0, 1].set_title('POLSPEVP2')
axs[1, 0].plot(POLSPEVP3[:, 0], POLSPEVP3[:, 8], 'r', label='castor')
axs[1, 0].plot(POLSPEVP3_dr[:, 0], POLSPEVP3_dr[:, 8], 'b--', label='castor-rfe')
axs[1, 0].legend(loc='lower right')
axs[1, 0].set_title('POLSPEVP3')
axs[1, 1].plot(POLSPEST[:, 0], POLSPEST[:, 8], 'r', label='castor')
axs[1, 1].plot(POLSPEST_dr[:, 0], POLSPEST_dr[:, 8], 'b--', label='castor-rfe')
axs[1, 1].legend(loc='lower right')
axs[1, 1].set_title('POLSPEST')

for ax in axs.flat:
    ax.set(xlabel='k', ylabel='fscore')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

#plt.show()
plt.savefig(current_dir + '/results/' + 'comparison5_castor_dr=0_1.png')
#############################################################################################





#############################################################################################
# getting best results of each algorithm
#kameris_POLSPEVP1_max_fscore = max(POLSPEVP1[:, 4])
#indices = np.where(POLSPEVP1[:, 4] == kameris_POLSPEVP1_max_fscore)
#kameris_min_k = min(indices[0])

#castor_POLSPEVP1_max_fscore = max(POLSPEVP1[:, 8])
#indices = np.where(POLSPEVP1[:, 4] == castor_POLSPEVP1_max_fscore)
#castor_min_k = min(indices[0])

castor_metrics = []
kameris_metrics = []
for i, dataset in enumerate(datasets):
    file_name = current_dir + '/results/' + dataset + '_dr=0.csv'
    data = np.genfromtxt(file_name, delimiter=',')

    kameris_max_fscore = max(data[:, 4])
    indices = np.where(data[:, 4] == kameris_max_fscore)
    kameris_min_k = min(indices[0])

    kameris_metrics.append([kameris_max_fscore, kameris_min_k])

    castor_max_fscore = max(data[:, 8])
    indices = np.where(data[:, 8] == castor_max_fscore)
    castor_min_k = min(indices[0])

    castor_metrics.append([castor_max_fscore, castor_min_k])

castor_metrics = np.array(castor_metrics)
kameris_metrics = np.array(kameris_metrics)

print(kameris_metrics[:, 0])
print(range(0, 8))

fig, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.plot(datasets, kameris_metrics[:, 0], 'r', label='kameris fscore')
ax1.legend(loc='lower right')
ax1.set(xlabel='datasets', ylabel='max fscore')

#ax1.set_title('POLSPEVP1')
ax2.plot(datasets, kameris_metrics[:, 1], 'r', label='kameris k-mer')
ax2.legend(loc='lower right')
ax2.set(xlabel='datasets', ylabel='k')
#ax2.set_xticklabels(['POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST', 'POLSPELT'])
#ax2.set_title('POLSPEVP1')
#axs[0, 1].plot(range(1, 9), kameris_metrics[:, 1], 'r', label='kameris k-mer')
#axs[0, 1].legend(loc='lower right')
#axs[0, 1].set(xlabel='datasets', ylabel='k')
#axs[0, 1].set_title('POLSPEVP1')

plt.xticks(rotation=45)
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

plt.savefig(current_dir + '/results/' + 'comparison_max_fscore.png')