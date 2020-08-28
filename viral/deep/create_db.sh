echo -e  "\n\nCreating EBOLA EBOSPECG..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/EBOLA/' '/home/vicente/datasets/MLDSP/' EBOSPECG kastor fasta

echo -e  "\n\nCreating HEPATITIS-B HBVGENCG..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/HEPATITIS-B/' '/home/vicente/datasets/MLDSP/' HBVGENCG kastor fasta


echo -e  "\n\nCreating HIV HIVGRPCG..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' '/home/vicente/datasets/MLDSP/' HIVGRPCG kastor fasta
echo -e  "\n\nCreating HIV HIVGRPCG..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' '/home/vicente/datasets/MLDSP/' HIVSUBCG kastor fasta
echo -e  "\n\nCreating HIV HIVGRPCG..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/HIV/' '/home/vicente/datasets/MLDSP/' HIVSUBPOL kastor fasta


echo -e  "\n\nCreating INFLUENZA INFSUBHA..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' '/home/vicente/datasets/MLDSP/' INFSUBHA kastor fasta
echo -e  "\n\nCreating INFLUENZA INFSUBMP..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' '/home/vicente/datasets/MLDSP/' INFSUBMP kastor fasta
echo -e  "\n\nCreating INFLUENZA INSUBFNA..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/INFLUENZA/' '/home/vicente/datasets/MLDSP/' INSUBFNA kastor fasta


echo -e  "\n\nCreating POLYOMAVIRUS POLSPELT..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/vicente/datasets/MLDSP/' POLSPELT kastor fasta
echo -e  "\n\nCreating POLYOMAVIRUS POLSPEST..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/vicente/datasets/MLDSP/' POLSPEST kastor fasta
echo -e  "\n\nCreating POLYOMAVIRUS POLSPEVP1..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/vicente/datasets/MLDSP/' POLSPEVP1 kastor fasta
echo -e  "\n\nCreating POLYOMAVIRUS POLSPEVP2..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/vicente/datasets/MLDSP/' POLSPEVP2 kastor fasta
echo -e  "\n\nCreating POLYOMAVIRUS POLSPEVP3..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/' '/home/vicente/datasets/MLDSP/' POLSPEVP3 kastor fasta


echo -e  "\n\nCreating PAPILLOMA HPVGENCG..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' '/home/vicente/datasets/MLDSP/' HPVGENCG kastor fasta
echo -e  "\n\nCreating PAPILLOMA HPVSPECG..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/' '/home/vicente/datasets/MLDSP/' HPVSPECG kastor fasta


echo -e  "\n\nCreating RHINOVIRUS RHISPECG..."
python3 buid_dataset_cgr.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/RHINOVIRUS/' '/home/vicente/datasets/MLDSP/' RHISPECG kastor fasta

echo -e  "\n\n FINISH DATASETS CREATION :) ###########################################################"
