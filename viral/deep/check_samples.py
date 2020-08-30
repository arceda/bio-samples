# este script cuenta la cantidad de muestras de las bases de datos

import sys
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO

path_dataset = "/home/siso/datasets/MLDSP/"
path_dataset = sys.argv[1]

datasets = glob.glob(path_dataset + "/*")

data_str = []

for db in datasets:
    print("\n" + db)

    train_labels = pd.read_csv(db + '/train_labels.csv', names=["sequence", "class"])  
    test_labels = pd.read_csv(db + '/test_labels.csv', names=["sequence", "class"])  

    seq_len_acc = 0
    for file in train_labels.values:
        sequences = SeqIO.parse( db + "/seq/" +  file[0], "fasta")
        for record in sequences:
            seq_len_acc +=  len(str(record.seq.upper()))

    for file in test_labels.values:
        sequences = SeqIO.parse( db + "/seq/" +  file[0], "fasta")
        for record in sequences:
            seq_len_acc +=  len(str(record.seq.upper()))

    avg_seq_len = seq_len_acc/ (train_labels.values.shape[0] + test_labels.values.shape[0])

    data_count_train = train_labels.groupby(['class']).count()
    data_count_test = test_labels.groupby(['class']).count()

    #print(data_count_train)
    #print(data_count_test)

    #print(data_count_train.values.shape[0], data_count_test.values.shape[0])

    if data_count_train.values.shape[0] != data_count_test.values.shape[0]:
        print("\n" + db, "have classes with 1 sample, impossible to train****************************")
        continue

    samples_per_class = data_count_train.values + data_count_test.values
    print("samples per class:", samples_per_class.T, "total samples:", np.sum(samples_per_class))

    file_name = db.split("/")[-1]
    data_str.append( [file_name, str(avg_seq_len), str(samples_per_class.shape[0]), str(samples_per_class.T), np.sum(samples_per_class)] )

print(data_str)

data_df = pd.DataFrame(data=data_str, columns=["db", "avg seq length", "num classes", "sample per class", "total"])  

print(data_df)

data_df.to_csv("db_description.csv")

# las bd que se estudiaron son:
# *Primates, Dengue, *Protist, Fungi, 
    

