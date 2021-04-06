#!/usr/bin/env python3
'''
A supervisory scripts to run all data through the
CUP primary tutmor type classificaiton pipeline.
Pipeline Steps:
    1. Ingest the data into the model.[]
    2. Preprocess the data using filtering.[X]
    3. Load the pretrained model.[X]
    4. Classify 32-way classification.[X]
    5. Provide the visualization of performances.[X]
    To do:
Portability:
	****Split the model building python script to a set of steps and break down into a snakeMake pipeline to execute each step
	This makes the model more dynamic and can include the raw or preprocessed data as an option for the clincian's labs*** []
	Build out the python script to make all the directories. [X]
	Build out the python script to run all of the options as flags. [X]
	Docker-SnakeMake-conda_env: model eval on local data set. [X]
	
Analysis:
	Build out the python script to run the model on the extended data. []
	
'''


import os
import re
import sys
import argparse
import subprocess


def parse_cli(args):
    '''
    Parses command line arguments.
    Parameters
    ----------
    args: list
        List of command line options followed by arguments.
    
    Returns
    -------
    args: dict
        Parsed arguments in key-value pairs.
    '''
    
    parser = argparse.ArgumentParser(
        description='A supervisory script to run all data through the'
                    + ' CUP classificaiton pipeline.'
    )
    required = parser.add_argument_group('required args')
    required.add_argument('-i', '--input', required=True,
                          help='Path to input directory containing'
                          		+ ' the methylation data.')
    required.add_argument('-o', '--output', required=True,
                          help='Path to output directory to save results.')
    required.add_argument('-m', '--model', required=True,
                          help='Path to CUP pretrained model.')
    required.add_argument('-c', '--classes_dir', required=True,
                          help='Path to classifications visualization data.')
    required.add_argument('-s', '--state', required=True,
                          help='Letting the model know if we have to analyze the data or just feature normalize')
    
    return vars(parser.parse_args(args))

def create_directories(output_dir_path):
    '''
    Creates directories to store pipeline results.
    Parameters
    ----------
    output_dir_path: string
        Path to output directory to make new directories.
    Returns
    -------
    paths: tuple of strings
        Paths of all the directories created. 
    '''

    model_path = os.path.join(output_dir_path, 'model_h5 files')
    classification_dir_path = os.path.join(output_dir_path, 
                                           'classification-results')

    os.makedirs(input_dir_path, exist_ok=True)
    os.makedirs(output_dir_path, exist_ok=True)
    os.makedirs(model_path, exist_ok=True)
    os.makedirs(classification_dir_path, exist_ok=True)
    
    paths = (
        model_path, 
        input_dir_path
        output_dir_path, 
        model_path,
        classification_dir_path
    )
    return paths
    

def main(cli_args=sys.argv[1:]):
    args = parse_cli(cli_args)
    input_dir_path = args['input']
    output_dir_path = args['output']
    model_path = args['model']
    classifications_dir_path = args['classes_dir']
    paths = create_directories(output_dir_path)
    
    input_dir_path = paths[0]
    output_dir_path = paths[1]
    model_path = paths[2]
    classification_dir_path = paths[3]
    state = args['state']

    for file_name in os.listdir(input_dir_path):
        if state == "Normalized":
		data_ingestion(input_dir_path)
	if state == "raw":
		analyze_methylation_data(input_dir_path)

    print('Feature Selection')
    feature_selection(important_features_dir_path, oversampled_dir_path,
                      feature_selected_dir_path)
    print('ANN processing your input')
    ANN(X)

if __name__ == '__main__':
    main()
