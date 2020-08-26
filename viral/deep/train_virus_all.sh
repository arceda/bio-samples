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


echo -e  "\n\nCreating INFLUENZA INFSUBHA..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' '/home/siso/datasets/MLDSP/' INFSUBHA kastor fasta
echo -e  "\n\nCreating INFLUENZA INFSUBMP..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' '/home/siso/datasets/MLDSP/' INFSUBMP kastor fasta
echo -e  "\n\nCreating INFLUENZA INSUBFNA..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' '/home/siso/datasets/MLDSP/' INSUBFNA kastor fasta


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


echo -e  "\n\nCreating RHINOVIRUS RHISPECG..."
python3 buid_dataset_cgr.py '/home/siso/projects/BIOINFORMATICS/datasets/VIRAL/RHINOVIRUS/' '/home/siso/datasets/MLDSP/' RHISPECG kastor fasta

echo -e  "\n\n FINISH DATASETS CREATION :) ###########################################################"




echo -e  "\n\ntraining EBOSPECG EBOSPECG dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' EBOSPECG  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' EBOSPECG  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' EBOSPECG

echo -e  "\n\ntraining HEPATITIS-B HBVGENCG dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HBVGENCG  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' HBVGENCG  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HBVGENCG

# HIV
echo -e  "\n\ntraining HIV HIVGRPCG dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HIVGRPCG  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVGRPCG  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HIVGRPCG

echo -e  "\n\ntraining HIV HIVSUBCG dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HIVSUBCG  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBCG  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HIVSUBCG

echo -e  "\n\ntraining HIV HIVSUBPOL dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HIVSUBPOL  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBPOL  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HIVSUBPOL



# POLIOMAVIRUS
echo -e  "\n\ntraining POLYOMAVIRUS POLSPELT dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' POLSPELT  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPELT  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' POLSPELT

echo -e  "\n\ntraining POLYOMAVIRUS POLSPEST dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' POLSPEST  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPEST  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' POLSPEST

echo -e  "\n\ntraining POLYOMAVIRUS POLSPEVP1 dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' POLSPEVP1  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPEVP1  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' POLSPEVP1

echo -e  "\n\ntraining POLYOMAVIRUS POLSPEVP2 dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' POLSPEVP2  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPEVP2  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' POLSPEVP2

echo -e  "\n\ntraining POLYOMAVIRUS  POLSPEVP3 dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' POLSPEVP3  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' POLSPEVP3  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' POLSPEVP3

# PAPILLOMA
echo -e  "\n\ntraining PAPILLOMA HPVGENCG dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HPVGENCG  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' HPVGENCG  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HPVGENCG

echo -e  "\n\ntraining PAPILLOMA HPVSPECG dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HPVSPECG  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' HPVSPECG  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HPVSPECG

# INFLUENZA
echo -e  "\n\ntraining INFLUENZA INFSUBHA dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' INFSUBHA  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBHA  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' INFSUBHA

echo -e  "\n\nraining INFLUENZA INFSUBMP dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' INFSUBMP  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBMP  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' INFSUBMP

echo -e  "\n\ntraining INFLUENZA INSUBFNA dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' INSUBFNA  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' INSUBFNA  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' INSUBFNA

# RHINOVIRUS
echo -e  "\n\ntraining RHINOVIRUS RHISPECG dataset..."
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' RHISPECG  0
python3 train_deep.py '/home/siso/datasets/MLDSP/' RHISPECG  50 128
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' RHISPECG
