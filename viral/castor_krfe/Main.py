#!/usr/bin/env python

###############
### Imports ###
############### 

from castor_krfe import data
from castor_krfe import extraction
from castor_krfe import evaluation
from castor_krfe import utils
from castor_krfe import __version__

import argparse

def parse_arguments():
    
    parser = argparse.ArgumentParser(description='CASTOR KRFE')

    # Modes
    parser.add_argument('-t', '--train',
                        action='store_true',
                        help="Extract best kmer features and"
                        "train a model with a dataset of labeled"
                        "genomic sequences")
 
    parser.add_argument('-e', '--evaluate',
                        action='store_true',
                        help="Validate dataset")
 
    parser.add_argument('-p', '--predict',
                        action='store_true',
                        help="Predict classes of genomic sequences")

    # Fasta file
    parser.add_argument('-f', '--fasta',
                        required=True,
                        help="Fasta file containing the sequences"
                        "The file is used in the three modes")

    # Class file
    parser.add_argument('-c', '--csv',
                        nargs='?',
                        help="CSV file containing the labels of the sequences"
                        "with the format : ID,label"
                        "The file could be used in the training/evaluation modes")

    # Model File
    parser.add_argument('-m', '--model',
                        nargs='?',
                        help="File of the model")
    
    # Extracted Kmer file
    parser.add_argument('-k', '--kmers',
                        nargs='?',
                        help="File of the extracted k-mers list")
    
    # Output folder
    parser.add_argument('-o', '--output',
                        nargs='?',
                        default='Output/',
                        help="Output directory. If it is not specified, "
                        "the program uses Output/ folder")
    
    # Parameters
    parser.add_argument('--threshold',
                        type=float,
                        nargs='?',
                        default=1,
                        help="Percentage of performance loss in terms of"
                        "F-measure to reduce the number of attributes")
 
    parser.add_argument('--kmin',
                        type=int,
                        nargs='?',
                        default=5,
                        help="Minimum length of k-mer(s)")

    parser.add_argument('--kmax',
                        type=int,
                        nargs='?',
                        default=5,
                        help="Maximum length of k-mer(s)")
    
    parser.add_argument('--fmin',
                        type=int,
                        nargs='?',
                        default=1,
                        help="Minimum number of features to identify")
    
    parser.add_argument('--fmax',
                        type=int,
                        nargs='?',
                        default=50,
                        help="Maximum number of features to identify")
    

    # Versoin
    parser.add_argument('--version',
                        action='version',
                        version=__version__)

    args = parser.parse_args()

    return args


def main():

    #######################
    ### Parse arguments ###
    #######################

    args = parse_arguments()

    # Modes
    mode_training = args.train
    mode_evaluation = args.evaluate
    mode_prediction = args.predict

    # Files and output directory
    kmers_file = args.kmers
    model_file = args.model

    class_file = args.csv
    fasta_file = args.fasta
    output_dir = args.output
    fig_file = output_dir + "Analysis"

    if not output_dir.endswith("/"): output_dir += "/"
    if not kmers_file: kmers_file = output_dir+'Kmers.txt'
    if not model_file: model_file = output_dir+'model.pkl'
    
    if mode_training or mode_evaluation:
        if not class_file:
            raise ValueError("CSV class file is required for this mode")

    # Parameters
    # Threshold (percentage of performance loss in terms of 
    # F-measure to reduce the number of attributes)
    T = args.threshold
    # Minimum length of k-mer(s)
    k_min = args.kmin
    # Maximum length of k-mer(s)
    k_max = args.kmax
    # Minimum number of features to identify
    features_min = args.fmin
    # Maximum number of features to identify
    features_max = args.fmax

    ####################
    ## Training mode ###
    ####################
    if mode_training:
        print("\nTraining mode\n")
        # Get training dataset
        print("Loading of the training dataset...")
        training_data = data.generateLabeledData(fasta_file, class_file)

        ###########################
        ### FEATURE EXTRACTION  ###
        ###########################
        best_k_mers, best_k_length = extraction.extractKmers(T, training_data,
                k_min, k_max, features_min, features_max, fig_file, kmers_file)
        print("Identified k =", best_k_length)
        print("Identified k-mers  =", best_k_mers)

        ##################
        ### EVALUATION ###
        ##################
        evaluation.cross_validation(training_data, best_k_mers)
        model = evaluation.train_model(training_data, best_k_mers) 
        utils.save_model(model, model_file)

    ######################
    ## Evaluation mode ###
    ######################
    if mode_evaluation:
        print("\nEvaluation mode\n")
        # Get testing dataset
        print("Loading of the testing dataset...")
        testing_data = data.generateLabeledData(fasta_file, class_file)
        
        # Load Model
        print("Loading the model from file...")
        model = utils.load_model(model_file)
        
        # Evaluation
        print("\nEvaluation")
        # the value of best_k_mers will be fetched from a file
        best_k_mers = utils.fetch_list_from_file(kmers_file)
        eval_file = output_dir+"Evaluation.txt"
        evaluation.evaluation(model, testing_data, best_k_mers, eval_file)
 
    ######################
    ## Prediction mode ###
    ######################
    if mode_prediction:
        print("\nPrediction mode\n")
        # Get testing dataset
        print("Loading of the testing dataset...")
        testing_data = data.generateData(fasta_file)
        
        # Load Model
        print("Loading the model from file...")
        model = utils.load_model(model_file)
        
        # Prediction
        print("\nPrediction")
        # the value of best_k_mers will be fetched from a file
        pred_file = output_dir+"Prediction.txt"
        best_k_mers = utils.fetch_list_from_file(kmers_file)
        evaluation.prediction(model, testing_data, best_k_mers, pred_file)


if __name__ == "__main__":
    main()
