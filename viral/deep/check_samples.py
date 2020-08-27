# este script cuenta la cantidad de muestras de las bases de datos

import sys
import glob


path_dataset = "'/home/siso/datasets/MLDSP/"
path_dataset = sys.argv[1]

datasets = glob.glob(path_dataset + "/*")

for db in datasets:
    print(db)