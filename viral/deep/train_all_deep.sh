#!/bin/bash
echo "welcome to vicente scripts"
echo "this script will train mldsp, cnn, kameris and castor over all datasets."

echo "training Protists dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Protists  0 50

echo "training Fungi dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Fungi  0 50

echo "training Plants dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Plants  0 50

echo "training Amphibians dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Amphibians  0 50

echo "training Mammals dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Mammals  0 50

echo "training Insects dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Insects  0 50

echo "training 3classes dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' 3classes  0 50

echo "training Vertebrates dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Vertebrates  0 50

