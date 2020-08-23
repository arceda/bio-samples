#!/bin/bash
echo "welcome to vicente scripts"
echo "this script will train mldsp, cnn, kameris and castor over all datasets."

echo "training Protists dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Protists  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Protists  
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Protists  

echo "training Fungi dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Fungi  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Fungi  
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Fungi  

echo "training Plants dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Plants  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Plants  
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Plants  

echo "training Amphibians dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Amphibians  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Amphibians  
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Amphibians  

echo "training Mammals dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Mammals  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Mammals  
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Mammals  

echo "training Insects dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Insects  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Insects  
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Insects  

echo "training 3classes dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' 3classes  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' 3classes  
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' 3classes  

echo "training Vertebrates dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Vertebrates  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' Vertebrates  
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Vertebrates  
