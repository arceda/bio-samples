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


#dataset_path = sys.argv[1]
#dataset_path = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/"


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

    #break


    
