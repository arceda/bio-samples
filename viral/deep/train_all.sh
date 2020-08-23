#!/bin/bash
echo "welcome to vicente scripts"
echo "this script will train mldsp, cnn, kameris and castor over all datasets."

echo "training Protists dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Protists  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Protists  0
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Protists  0

echo "training Fungi dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Fungi  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Fungi  0
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Fungi  0

echo "training Plants dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Plants  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Plants  0
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Plants  0

echo "training Amphibians dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Amphibians  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Amphibians  0
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Amphibians  0

echo "training Mammals dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Mammals  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Mammals  0
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Mammals  0

echo "training Insects dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Insects  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Insects  0
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Insects  0

echo "training 3classes dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' 3classes  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' 3classes  0
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' 3classes  0

echo "training Vertebrates dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Vertebrates  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Vertebrates  0
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Vertebrates  0
