#!/bin/bash
echo -e  "welcome to siso scripts"
echo -e  "this script will train mldsp, cnn, kameris and castor over all datasets."

# database creation
echo -e  "\n\nCreating datasets..."

echo -e  "\n\nCreating EBOLA EBOSPECG..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/EBOLA/' '/home/siso/datasets/MLDSP/' EBOSPECG kastor fasta


echo -e  "\n\nCreating HEPATITIS-B HBVGENCG..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HEPATITIS-B/' '/home/siso/datasets/MLDSP/' HBVGENCG kastor fasta


echo -e  "\n\nCreating HIV HIVGRPCG..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' '/home/siso/datasets/MLDSP/' HIVGRPCG kastor fasta
echo -e  "\n\nCreating HIV HIVGRPCG..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' '/home/siso/datasets/MLDSP/' HIVSUBCG kastor fasta
echo -e  "\n\nCreating HIV HIVGRPCG..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' '/home/siso/datasets/MLDSP/' HIVSUBPOL kastor fasta


#echo -e  "\n\nCreating INFLUENZA INFSUBHA..."
#python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' '/home/siso/datasets/MLDSP/' INFSUBHA kastor fasta
#echo -e  "\n\nCreating INFLUENZA INFSUBMP..."
#python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' '/home/siso/datasets/MLDSP/' INFSUBMP kastor fasta
#echo -e  "\n\nCreating INFLUENZA INSUBFNA..."
#python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' '/home/siso/datasets/MLDSP/' INSUBFNA kastor fasta


echo -e  "\n\nCreating POLYOMAVIRUS POLSPELT..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/siso/datasets/MLDSP/' POLSPELT kastor fasta
echo -e  "\n\nCreating POLYOMAVIRUS POLSPEST..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/siso/datasets/MLDSP/' POLSPEST kastor fasta
echo -e  "\n\nCreating POLYOMAVIRUS POLSPEVP1..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/siso/datasets/MLDSP/' POLSPEVP1 kastor fasta
echo -e  "\n\nCreating POLYOMAVIRUS POLSPEVP2..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/siso/datasets/MLDSP/' POLSPEVP2 kastor fasta
echo -e  "\n\nCreating POLYOMAVIRUS POLSPEVP3..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/siso/datasets/MLDSP/' POLSPEVP3 kastor fasta


echo -e  "\n\nCreating PAPILLOMA HPVGENCG..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' '/home/siso/datasets/MLDSP/' HPVGENCG kastor fasta
echo -e  "\n\nCreating PAPILLOMA HPVSPECG..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' '/home/siso/datasets/MLDSP/' HPVSPECG kastor fasta

echo -e  "\n\n FINISH DATASETS CREATION :) ###########################################################"






echo -e  "\n\ntraining EBOSPECG EBOSPECG dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/EBOLA/' EBOSPECG  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/EBOLA/' EBOSPECG  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/EBOLA/' EBOSPECG  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/EBOLA/' EBOSPECG

echo -e  "\n\ntraining HEPATITIS-B HBVGENCG dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HEPATITIS-B/' HBVGENCG  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HEPATITIS-B/' HBVGENCG  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HEPATITIS-B/' HBVGENCG  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HEPATITIS-B/' HBVGENCG

# HIV
echo -e  "\n\ntraining HIV HIVGRPCG dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVGRPCG  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVGRPCG  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVGRPCG  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVGRPCG

echo -e  "\n\ntraining HIV HIVSUBCG dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVSUBCG  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVSUBCG  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVSUBCG  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVSUBCG

echo -e  "\n\ntraining HIV HIVSUBPOL dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVSUBPOL  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVSUBPOL  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVSUBPOL  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' HIVSUBPOL



# POLIOMAVIRUS
echo -e  "\n\ntraining POLYOMAVIRUS POLSPELT dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPELT  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPELT  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPELT  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPELT

echo -e  "\n\ntraining POLYOMAVIRUS POLSPEST dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEST  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEST  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEST  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEST

echo -e  "\n\ntraining POLYOMAVIRUS POLSPEVP1 dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP1  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP1  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP1  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP1

echo -e  "\n\ntraining POLYOMAVIRUS POLSPEVP2 dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP2  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP2  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP2  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP2

echo -e  "\n\ntraining POLYOMAVIRUS  POLSPEVP3 dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP3  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP3  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP3  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' POLSPEVP3

# PAPILLOMA
echo -e  "\n\ntraining PAPILLOMA HPVGENCG dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' HPVGENCG  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' HPVGENCG  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' HPVGENCG  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' HPVGENCG

echo -e  "\n\ntraining PAPILLOMA HPVSPECG dataset..."
python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' HPVSPECG  0
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' HPVSPECG  10 128
python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' HPVSPECG  50 128
python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' HPVSPECG

# INFLUENZA
#echo -e  "\n\ntraining INFLUENZA INFSUBHA dataset..."
#python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INFSUBHA  0
#python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INFSUBHA  10 128
#python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INFSUBHA  50 128
#python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INFSUBHA

#echo -e  "\n\nraining INFLUENZA INFSUBMP dataset..."
#python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INFSUBMP  0
#python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INFSUBMP  10 128
#python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INFSUBMP  50 128
#python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INFSUBMP

#echo -e  "\n\ntraining INFLUENZA INSUBFNA dataset..."
#python3 train_mldsp.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INSUBFNA  0
#python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INSUBFNA  10 128
#python3 train_deep.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INSUBFNA  50 128
#python3 train_kameris_castor.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' INSUBFNA

