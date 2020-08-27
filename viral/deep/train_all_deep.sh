#!/bin/bash
echo -e "welcome to vicente scripts"
echo -e "this script will train mldsp, cnn, kameris and castor over all datasets."

echo -e "training Primates dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Primates  50 128

echo -e "training Dengue dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Dengue  50 128

echo -e "training Protists dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Protists  50 128

echo -e "training Fungi dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Fungi  50 128

echo -e "training Plants dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Plants  50 128

echo -e "training Amphibians dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Amphibians  50 128

echo -e "training Mammals dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Mammals  50 128

echo -e "training Insects dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Insects  50 128

echo -e "training 3classes dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' 3classes  50 128

echo -e "training Vertebrates dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Vertebrates  50 128




echo -e  "\n\ntraining EBOSPECG EBOSPECG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' EBOSPECG  50 128

echo -e  "\n\ntraining HEPATITIS-B HBVGENCG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HBVGENCG  50 128



echo -e  "\n\ntraining HIV HIVGRPCG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVGRPCG  50 128

echo -e  "\n\ntraining HIV HIVSUBCG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBCG  50 128

echo -e  "\n\ntraining HIV HIVSUBPOL dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBPOL  50 128


# POLIOMAVIRUS
echo -e  "\n\ntraining POLYOMAVIRUS POLSPELT dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPELT  50 128

echo -e  "\n\ntraining POLYOMAVIRUS POLSPEST dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPEST  50 128

echo -e  "\n\ntraining POLYOMAVIRUS POLSPEVP1 dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPEVP1  50 128

echo -e  "\n\ntraining POLYOMAVIRUS POLSPEVP2 dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPEVP2  50 128

echo -e  "\n\ntraining POLYOMAVIRUS  POLSPEVP3 dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPEVP3  50 128

# PAPILLOMA
echo -e  "\n\ntraining PAPILLOMA HPVGENCG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HPVGENCG  50 128

echo -e  "\n\ntraining PAPILLOMA HPVSPECG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HPVSPECG  50 128

# INFLUENZA
echo -e  "\n\ntraining INFLUENZA INFSUBHA dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBHA  50 128

echo -e  "\n\nraining INFLUENZA INFSUBMP dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBMP  50 128

echo -e  "\n\ntraining INFLUENZA INSUBFNA dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' INSUBFNA  50 128

# RHINOVIRUS
echo -e  "\n\ntraining RHINOVIRUS RHISPECG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' RHISPECG  50 128
