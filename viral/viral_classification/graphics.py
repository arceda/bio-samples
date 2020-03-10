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
#datasets = ['HIVGRPCG', 'HIVSUBCG', 'HIVSUBPOL'] 

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

'''
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
HIVGRPCG_dr = np.genfromtxt(current_dir + '/results/HIVGRPCG_dr=1.csv', delimiter=',')
HIVSUBCG_dr = np.genfromtxt(current_dir + '/results/HIVSUBCG_dr=1.csv', delimiter=',')
HIVSUBPOL_dr = np.genfromtxt(current_dir + '/results/HIVSUBPOL_dr=1.csv', delimiter=',')

fscore_kameris_pos = 4
fscore_castor_pos = 8
k_pos = 0
'''

POLSPEVP1 = np.genfromtxt(current_dir + '/results/POLSPEVP1_dr=0_nfeatures.csv', delimiter=',')
POLSPEVP2 = np.genfromtxt(current_dir + '/results/POLSPEVP2_dr=0_nfeatures.csv', delimiter=',')
POLSPEVP3 = np.genfromtxt(current_dir + '/results/POLSPEVP3_dr=0_nfeatures.csv', delimiter=',')
POLSPEST = np.genfromtxt(current_dir + '/results/POLSPEST_dr=0_nfeatures.csv', delimiter=',')

POLSPELT = np.genfromtxt(current_dir + '/results/POLSPELT_dr=0_nfeatures.csv', delimiter=',')
HIVGRPCG = np.genfromtxt(current_dir + '/results/HIVGRPCG_dr=0_nfeatures.csv', delimiter=',')
HIVSUBCG = np.genfromtxt(current_dir + '/results/HIVSUBCG_dr=0_nfeatures.csv', delimiter=',')
HIVSUBPOL = np.genfromtxt(current_dir + '/results/HIVSUBPOL_dr=0_nfeatures.csv', delimiter=',')

POLSPEVP1_dr = np.genfromtxt(current_dir + '/results/POLSPEVP1_dr=1_nfeatures.csv', delimiter=',')
POLSPEVP2_dr = np.genfromtxt(current_dir + '/results/POLSPEVP2_dr=1_nfeatures.csv', delimiter=',')
POLSPEVP3_dr = np.genfromtxt(current_dir + '/results/POLSPEVP3_dr=1_nfeatures.csv', delimiter=',')
POLSPEST_dr = np.genfromtxt(current_dir + '/results/POLSPEST_dr=1_nfeatures.csv', delimiter=',')

POLSPELT_dr = np.genfromtxt(current_dir + '/results/POLSPELT_dr=1_nfeatures.csv', delimiter=',')
HIVGRPCG_dr = np.genfromtxt(current_dir + '/results/HIVGRPCG_dr=1_nfeatures.csv', delimiter=',')
HIVSUBCG_dr = np.genfromtxt(current_dir + '/results/HIVSUBCG_dr=1_nfeatures.csv', delimiter=',')
HIVSUBPOL_dr = np.genfromtxt(current_dir + '/results/HIVSUBPOL_dr=1_nfeatures.csv', delimiter=',')

fscore_kameris_pos = 4
fscore_castor_pos = 9
k_pos = 0

############################################################################################
# comparison between 'POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST' without dimentionality reduction
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(POLSPEVP1[:, k_pos], POLSPEVP1[:, fscore_kameris_pos], 'r.-', label='kameris')
axs[0, 0].plot(POLSPEVP1[:, k_pos], POLSPEVP1[:, fscore_castor_pos], 'b', label='castor')
axs[0, 0].plot(POLSPEVP1_dr[:, k_pos], POLSPEVP1_dr[:, fscore_kameris_pos], 'g--', label='kameris-sdv')
axs[0, 0].plot(POLSPEVP1_dr[:, k_pos], POLSPEVP1_dr[:, fscore_castor_pos], 'm-.', label='castor-krfe')
axs[0, 0].legend(loc='lower right', fontsize=8)
axs[0, 0].set_title('POLSPEVP1')
axs[0, 1].plot(POLSPEVP2[:, k_pos], POLSPEVP2[:, fscore_kameris_pos], 'r.-', label='kameris')
axs[0, 1].plot(POLSPEVP2[:, k_pos], POLSPEVP2[:, fscore_castor_pos], 'b', label='castor')
axs[0, 1].plot(POLSPEVP2_dr[:, k_pos], POLSPEVP2_dr[:, fscore_kameris_pos], 'g--', label='kameris-sdv')
axs[0, 1].plot(POLSPEVP2_dr[:, k_pos], POLSPEVP2_dr[:, fscore_castor_pos], 'm-.', label='castor-krfe')
axs[0, 1].legend(loc='lower right', fontsize=8)
axs[0, 1].set_title('POLSPEVP2')
axs[1, 0].plot(POLSPEVP3[:, k_pos], POLSPEVP3[:, fscore_kameris_pos], 'r.-', label='kameris')
axs[1, 0].plot(POLSPEVP3[:, k_pos], POLSPEVP3[:, fscore_castor_pos], 'b', label='castor')
axs[1, 0].plot(POLSPEVP3_dr[:, k_pos], POLSPEVP3_dr[:, fscore_kameris_pos], 'g--', label='kameris-sdv')
axs[1, 0].plot(POLSPEVP3_dr[:, k_pos], POLSPEVP3_dr[:, fscore_castor_pos], 'm-.', label='castor-krfe')
axs[1, 0].legend(loc='lower right', fontsize=8)
axs[1, 0].set_title('POLSPEVP3')
axs[1, 1].plot(POLSPEST[:, k_pos], POLSPEST[:, fscore_kameris_pos], 'r.-', label='kameris')
axs[1, 1].plot(POLSPEST[:, k_pos], POLSPEST[:, fscore_castor_pos], 'b', label='castor')
axs[1, 1].plot(POLSPEST_dr[:, k_pos], POLSPEST_dr[:, fscore_kameris_pos], 'g--', label='kameris-sdv')
axs[1, 1].plot(POLSPEST_dr[:, k_pos], POLSPEST_dr[:, fscore_castor_pos], 'm-.', label='castor-krfe')
axs[1, 1].legend(loc='lower right', fontsize=8)
axs[1, 1].set_title('POLSPEST')

for ax in axs.flat:
    ax.set(xlabel='k', ylabel='fscore')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

#plt.show()
plt.savefig(current_dir + '/results/' + 'comparison3_kameris_castor_dr=0_1.png', dpi = 300)
#############################################################################################


############################################################################################
# comparison between 'POLSPELT', 'HIVGRPCG', 'HIVSUBCG', 'HIVSUBPOL' without dimentionality reduction
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(POLSPELT[:, k_pos], POLSPELT[:, fscore_kameris_pos], 'r.-', label='kameris')
axs[0, 0].plot(POLSPELT[:, k_pos], POLSPELT[:, fscore_castor_pos], 'b', label='castor')
axs[0, 0].plot(POLSPELT_dr[:, k_pos], POLSPELT_dr[:, fscore_kameris_pos], 'g--', label='kameris-sdv')
axs[0, 0].plot(POLSPELT_dr[:, k_pos], POLSPELT_dr[:, fscore_castor_pos], 'm-.', label='castor-krfe')
axs[0, 0].legend(loc='lower right', fontsize=8)
axs[0, 0].set_title('POLSPELT')
axs[0, 1].plot(HIVGRPCG[:, k_pos], HIVGRPCG[:, fscore_kameris_pos], 'r.-', label='kameris')
axs[0, 1].plot(HIVGRPCG[:, k_pos], HIVGRPCG[:, fscore_castor_pos], 'b', label='castor')
axs[0, 1].plot(HIVGRPCG_dr[:, k_pos], HIVGRPCG_dr[:, fscore_kameris_pos], 'g--', label='kameris-sdv')
axs[0, 1].plot(HIVGRPCG_dr[:, k_pos], HIVGRPCG_dr[:, fscore_castor_pos], 'm-.', label='castor-krfe')
axs[0, 1].legend(loc='lower right', fontsize=8)
axs[0, 1].set_title('HIVGRPCG')
axs[1, 0].plot(HIVSUBCG[:, k_pos], HIVSUBCG[:, fscore_kameris_pos], 'r.-', label='kameris')
axs[1, 0].plot(HIVSUBCG[:, k_pos], HIVSUBCG[:, fscore_castor_pos], 'b', label='castor')
axs[1, 0].plot(HIVSUBCG_dr[:, k_pos], HIVSUBCG_dr[:, fscore_kameris_pos], 'g--', label='kameris-sdv')
axs[1, 0].plot(HIVSUBCG_dr[:, k_pos], HIVSUBCG_dr[:, fscore_castor_pos], 'm-.', label='castor-krfe')
axs[1, 0].legend(loc='lower right', fontsize=8)
axs[1, 0].set_title('HIVSUBCG')
axs[1, 1].plot(HIVSUBPOL[:, k_pos], HIVSUBPOL[:, fscore_kameris_pos], 'r.-', label='kameris')
axs[1, 1].plot(HIVSUBPOL[:, k_pos], HIVSUBPOL[:, fscore_castor_pos], 'b', label='castor')
axs[1, 1].plot(HIVSUBPOL_dr[:, k_pos], HIVSUBPOL_dr[:, fscore_kameris_pos], 'g--', label='kameris-sdv')
axs[1, 1].plot(HIVSUBPOL_dr[:, k_pos], HIVSUBPOL_dr[:, fscore_castor_pos], 'm-.', label='castor-krfe')
axs[1, 1].legend(loc='lower right', fontsize=8)
axs[1, 1].set_title('HIVSUBPOL')

for ax in axs.flat:
    ax.set(xlabel='k', ylabel='fscore')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

#plt.show()
plt.savefig(current_dir + '/results/' + 'comparison4_kameris_castor_dr=0_1.png', dpi = 300)
#############################################################################################

'''
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

'''



#############################################################################################
# getting best results of each algorithm
#kameris_POLSPEVP1_max_fscore = max(POLSPEVP1[:, 4])
#indices = np.where(POLSPEVP1[:, 4] == kameris_POLSPEVP1_max_fscore)
#kameris_min_k = min(indices[0])

#castor_POLSPEVP1_max_fscore = max(POLSPEVP1[:, 8])
#indices = np.where(POLSPEVP1[:, 4] == castor_POLSPEVP1_max_fscore)
#castor_min_k = min(indices[0])

castor_metrics = []
castor_metrics_dr = []
kameris_metrics = []
kameris_metrics_dr = []
for i, dataset in enumerate(datasets):
    
    file_name = current_dir + '/results/' + dataset + '_dr=0_nfeatures.csv'
    data = np.genfromtxt(file_name, delimiter=',')

    kameris_max_fscore = max(data[:, 4])
    indices = np.where(data[:, 4] == kameris_max_fscore)
    kameris_min_k_index = min(indices[0])
    kameris_min_k = data[kameris_min_k_index, 0] # in this column is k value

    kameris_metrics.append([kameris_max_fscore, kameris_min_k])

    castor_max_fscore = max(data[:, 9])
    indices = np.where(data[:, 9] == castor_max_fscore)
    castor_min_k_index = min(indices[0])
    castor_min_k = data[castor_min_k_index, 0] # in this column is k value

    castor_metrics.append([castor_max_fscore, castor_min_k])



    file_name = current_dir + '/results/' + dataset + '_dr=1_nfeatures.csv'
    data = np.genfromtxt(file_name, delimiter=',')

    kameris_max_fscore = max(data[:, 4])
    indices = np.where(data[:, 4] == kameris_max_fscore)
    kameris_min_k_index = min(indices[0])
    kameris_min_k = data[kameris_min_k_index, 0] # in this column is k value

    kameris_metrics_dr.append([kameris_max_fscore, kameris_min_k])

    castor_max_fscore = max(data[:, 9])
    indices = np.where(data[:, 9] == castor_max_fscore)
    castor_min_k_index = min(indices[0])
    castor_min_k = data[castor_min_k_index, 0] # in this column is k value

    castor_metrics_dr.append([castor_max_fscore, castor_min_k])

castor_metrics = np.array(castor_metrics)
kameris_metrics = np.array(kameris_metrics)
castor_metrics_dr = np.array(castor_metrics_dr)
kameris_metrics_dr = np.array(kameris_metrics_dr)



fig, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.plot(datasets, kameris_metrics[:, 0], 'r.-', label='kameris')
ax1.plot(datasets, castor_metrics[:, 0], 'b', label='castor')
ax1.plot(datasets, kameris_metrics_dr[:, 0], 'g--', label='kameris-svd')
ax1.plot(datasets, castor_metrics_dr[:, 0], 'm-.', label='castor-rfe')
ax1.legend(loc='lower left', fontsize=8)
ax1.set(xlabel='datasets', ylabel='max fscore')
ax1.label_outer()

#ax1.set_title('POLSPEVP1')
ax2.plot(datasets, kameris_metrics[:, 1], 'r.-', label='kameris')
ax2.plot(datasets, castor_metrics[:, 1], 'b', label='castor')
ax2.plot(datasets, kameris_metrics_dr[:, 1], 'g--', label='kameris-svd')
ax2.plot(datasets, castor_metrics_dr[:, 1], 'm-.', label='castor-rfe')
ax2.legend(loc='upper left', fontsize=8)
ax2.set(xlabel='datasets', ylabel='k')
ax1.label_outer()
#ax2.set_xticklabels(['POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST', 'POLSPELT'])
#ax2.set_title('POLSPEVP1')
#axs[0, 1].plot(range(1, 9), kameris_metrics[:, 1], 'r', label='kameris k-mer')
#axs[0, 1].legend(loc='lower right')
#axs[0, 1].set(xlabel='datasets', ylabel='k')
#axs[0, 1].set_title('POLSPEVP1')

plt.xticks(rotation=45, fontsize=6)
# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()

plt.savefig(current_dir + '/results/' + 'comparison_max_fscore.png')











#################################################################################################
############## bar ##############################################################################
#################################################################################################
'''
datasets = ['POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST', 'POLSPELT'] 
kameris_nfeatures_list = []
castor_nfeatures_list = []

for i, dataset in enumerate(datasets):
    file_name = current_dir + '/results/' + dataset + '_dr=0_nfeatures.csv'
    data = np.genfromtxt(file_name, delimiter=',')

    kameris_nfeatures = data[:, 5]
    castor_nfeatures = data[:, 10]

    kameris_nfeatures_list.append(kameris_nfeatures)
    castor_nfeatures_list.append(castor_nfeatures)

    

x = np.arange(len(datasets))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, kameris_nfeatures_list, width, label='kameris')
rects2 = ax.bar(x + width/2, castor_nfeatures_list, width, label='castor')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('nfeatures')
#ax.set_title('Scores by group and gender')
ax.set_xticks(x)
ax.set_xticklabels(datasets)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)

fig.tight_layout()

plt.savefig(current_dir + '/results/' + 'comparison_nfeatures.png')

'''








