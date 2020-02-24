import numpy as np
import sys
import csv
import matplotlib as mpl 
## agg backend is used to create plot as a .png file
#mpl.use('agg')
import matplotlib.pyplot as plt 
from matplotlib import pyplot

#datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 
#            'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT', 'HEPATITIS-B/HBVGENCG'] 
datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 
            'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT'] 
#datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2'] 

dataset_path = sys.argv[1]
#dataset_path = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/"


def plot_box_plot(data, ylabel):
    # Create a figure instance
    fig = plt.figure(1, figsize=(9, 6))

    # Create an axes instance
    ax = fig.add_subplot(111)

    #ax.set_title('Accuracy')
    ax.set_axisbelow(True)
    ax.set_xlabel('Datasets')
    ax.set_ylabel(ylabel)

    # Create the boxplot
    bp = ax.boxplot(data, showmeans=True)

    #ax.set_ylim(0, 1.2)
    #ax.set_xticklabels(['POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST', 'POLSPELT', 'HBVGENCG'], rotation=45, fontsize=8)
    ax.set_xticklabels(['POLSPEVP1', 'POLSPEVP2', 'POLSPEVP3', 'POLSPEST', 'POLSPELT'])
    # Save the figure
    #fig.savefig('viral/viral_classification/results/' + ylabel + '.png', bbox_inches='tight')
    plt.show()
    ax.clear()

total_acc = []
total_precision = []
total_recall = []
total_fscore = []
total_best_k = []
total_best_nfeatures = []

for i, dataset in enumerate(datasets):
    print(i, "\n\nEVALUATING DATASET: ", dataset, "...")
    acc = np.array([])
    precision = np.array([])
    recall = np.array([])
    fscore = np.array([])
    best_k = np.array([])
    best_nfeatures = np.array([])
    k_mers = np.array([])
    times = np.array([])



    #results = np.genfromtxt(dataset_path + dataset + '/results.metrics', delimiter=';')
    with open(dataset_path + dataset + '/results.metrics', 'r') as file:
        reader = csv.reader(file, delimiter = ';')
        next(reader) #skip the first row
        for row in reader:           
            db_name = row[0]
            acc             = np.append(acc, float(row[1]))
            precision       = np.append(precision, float(row[2]))
            recall          = np.append(recall, float(row[3]))
            fscore          = np.append(fscore, float(row[4]))
            best_k          = np.append(best_k, int(row[5]))
            best_nfeatures  = np.append(best_nfeatures, int(row[6]))

    pyplot.hist(best_nfeatures)
    pyplot.show()


    total_acc.append(acc)   
    total_precision.append(precision)
    total_recall.append(recall)
    total_fscore.append(fscore)
    total_best_k.append(best_k)
    total_best_nfeatures.append(best_nfeatures)

plot_box_plot(total_acc, 'acc')
plot_box_plot(total_precision, 'precision')
plot_box_plot(total_recall, 'recall')
plot_box_plot(total_fscore, 'fscore')
plot_box_plot(total_best_k, 'k')
plot_box_plot(total_best_nfeatures, 'nfeatures')

