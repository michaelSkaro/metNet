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
	Build out the python script to run the model on the extended data. [X]
	
'''


import click
import io
import math
import os
import re
import time
import warnings
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost
from matplotlib import pyplot
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import explained_variance_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from xgboost import XGBRegressor
import collections
import random
import tensorflow as tf
print("Num GPUs Available: ", len(tf.config.list_physical_devices("GPU")))
import imblearn
from sklearn.model_selection import cross_val_score
from tensorflow import keras
from tensorflow.keras import Sequential, backend
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
import sklearn.metrics as metrics
warnings.filterwarnings("ignore")




@click.command()
@click.option('--path', default=1, help='PLease provide the path for the training data.')

# bunch of the ugly stuff
def data_preprocessing(X):
    # if the file has rowIDs, we can drop these
    if "Composite Element REF" in X.columns:
        X = X.drop("Composite Element REF", axis=1)
    if "Unnamed: 0" in X.columns:
        X = X.drop("Unnamed: 0", axis=1)
    # subset out the rows that have TCGA-OV
    X = X[X.cancerType != "TCGA-OV"]
    # X = X[X.cancerType != "TCGA-CHOL"]
    # combined READ and COAD
    X = X.replace({"cancerType": "TCGA-READ"}, "TCGA-COADREAD")
    X = X.replace({"cancerType": "TCGA-COAD"}, "TCGA-COADREAD")
    # substring the TCGA- out of the column

    # get rid of TCGA- in cancer type
    X["cancerType"] = X["cancerType"].str[5:]
    y = X["cancerType"]
    X = X.drop("cancerType", axis=1)
    X = pd.get_dummies(X, columns=["labels"])
    return X, y

def encode_y(y):
    y = np.array(y)
    y = y.reshape(-1, 1)
    encoder = OneHotEncoder(sparse=False)
    encoder.fit(y)
    y_encoded = encoder.transform(y)

    return y_encoded
    
    
seed = 19  # always 19 love you bubba
np.random.seed(seed)
## Artificial Neural Network
def build_classifier():
    classifier = Sequential()  # use the sequential layer
    classifier.add(
        Dense(units=80, kernel_initializer="uniform", activation="tanh", input_dim=4767)
    )
    classifier.add(Dropout(rate=0.5))
    classifier.add(
        Dense(units=80, kernel_initializer="uniform", activation="tanh", input_dim=4767)
    )
    classifier.add(Dropout(rate=0.5))
    classifier.add(
        Dense(units=80, kernel_initializer="uniform", activation="tanh", input_dim=4767)
    )
    ## here is the output layer, use softtmax for k classes > 2. Use signmoid for binary.
    classifier.add(Dense(units=30, kernel_initializer="uniform", activation="softmax"))
    opt = keras.optimizers.SGD(learning_rate=0.001)
    # Compiling the ANN
    classifier.compile(
        optimizer=opt, loss="categorical_crossentropy", metrics=["accuracy"]
    )
    
    return classifier

    
    

if __name__ == '__main__':
	X = pd.read_csv("Selected_features_all_methylation_data_TCGA.csv")
	X, y = data_preprocessing(X)
	y_encoded = encode_y(y)
	X_train, X_test, y_train, y_test = train_test_split(
    X, y_encoded, test_size=0.2, random_state=19, shuffle=True)  # always 19 <3 u bubba
    clf_nn = build_classifier()
    print(rounded_labels.shape)
    clf_nn.fit(X_train, y_train, batch_size=50, epochs=1200)
    tf.keras.models.save_model(clf_nn, filepath = 'myModel.h5')

	y_pred_nn = clf_nn.predict(X_test)
	y_pred = np.argmax(y_pred_nn, axis=1)
	rounded_labels = np.argmax(y_test, axis=1)
	print(rounded_labels)
    
    
    
    
    
    
    
    
