echo -e "welcome to vicente scripts"
echo -e "this script will train mldsp, cnn, kameris and castor over all datasets."

echo -e "training Primates dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Primates  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' Primates  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' Primates  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Primates 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Primates  

echo -e "training Dengue dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Dengue  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' Dengue  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' Dengue  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Dengue 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Dengue 

echo -e "training Protists dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Protists  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' Protists  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' Protists  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Protists 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Protists 

echo -e "training Fungi dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Fungi  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' Fungi  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' Fungi  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Fungi 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Fungi 

echo -e "training Plants dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Plants  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' Plants  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' Plants  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Plants 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Plants 

echo -e "training Amphibians dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Amphibians  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' Amphibians  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' Amphibians  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Amphibians 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Amphibians 

echo -e "training Insects dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Insects  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' Insects  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' Insects  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Insects 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Insects 

echo -e "training 3classes dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' 3classes  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' 3classes  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' 3classes  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' 3classes 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' 3classes 

echo -e "training Vertebrates dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' Vertebrates  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' Vertebrates  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' Vertebrates  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' Vertebrates 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' Vertebrates 


###################################################################################################
###################################################################################################






echo -e  "\n\ntraining HIV HIVGRPCG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVGRPCG  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVGRPCG  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVGRPCG  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HIVGRPCG 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HIVGRPCG 

echo -e  "\n\ntraining HIV HIVSUBCG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBCG  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBCG  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBCG  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HIVSUBCG 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HIVSUBCG 

echo -e  "\n\ntraining HIV HIVSUBPOL dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBPOL  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBPOL  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' HIVSUBPOL  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HIVSUBPOL 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HIVSUBPOL 





# INFLUENZA
echo -e  "\n\ntraining INFLUENZA INFSUBHA dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBHA  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBHA  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBHA  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' INFSUBHA 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' INFSUBHA 

echo -e  "\n\nraining INFLUENZA INFSUBMP dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBMP  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBMP  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' INFSUBMP  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' INFSUBMP 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' INFSUBMP 

echo -e  "\n\ntraining INFLUENZA INSUBFNA dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' INSUBFNA  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' INSUBFNA  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' INSUBFNA  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' INSUBFNA 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' INSUBFNA 

# EBOLA
echo -e  "\n\ntraining EBOLA EBOSPECG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' EBOSPECG  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' EBOSPECG  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' EBOSPECG  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' EBOSPECG 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' EBOSPECG 

# HEPATITIS-B
echo -e  "\n\ntraining HEPATITIS-B HBVGENCG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HBVGENCG  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' HBVGENCG  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' HBVGENCG  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HBVGENCG 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HBVGENCG 


# RHINOVIRUS
echo -e  "\n\ntraining RHINOVIRUS RHISPECG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' RHISPECG  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' RHISPECG  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' RHISPECG  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' RHISPECG 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' RHISPECG 

# PAPILLOMA
echo -e  "\n\ntraining PAPILLOMA HPVGENCG dataset..."
python3 train_deep.py '/home/siso/datasets/MLDSP/' HPVGENCG  100 128 tiny
python3 train_deep.py '/home/siso/datasets/MLDSP/' HPVGENCG  100 128 medium
python3 train_deep.py '/home/siso/datasets/MLDSP/' HPVGENCG  100 128 complex
python3 train_mldsp.py '/home/siso/datasets/MLDSP/' HPVGENCG 
python3 train_kameris_castor.py '/home/siso/datasets/MLDSP/' HPVGENCG 
