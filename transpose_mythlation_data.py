# import libraries for the gradient boosting
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xgboost
import math
from sklearn.model_selection import train_test_split
from sklearn.metrics import explained_variance_score
from xgboost import XGBRegressor
from matplotlib import pyplot
import re
import io
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import OneHotEncoder
import warnings
warnings.filterwarnings("ignore")

''' Make the reader for feature subsetting

Read in each dataframe, subset into the current row block 5000.
transform the DF to have features across the columns with all patients as row
Add a label to the edge. save file. 
'''

path = "/home/jovyan/storage/CUP/tsv/Compressed/*.gz"


def chunkIt(seq, num):
    '''
    Gonna split the feature sets into chunks
    '''
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def read_cut(file, start, end):
    sampleIDs = pd.read_table(file, delimiter="\t", compression="gzip", nrows=0).columns
    dat = pd.read_table(file, delimiter="\t", compression="gzip", skiprows=start, nrows=end - start, names=sampleIDs).fillna(0).set_index("Composite Element REF").T
    
    return dat

def CT(file):
    '''
    Str extract the Cancer type for the 
    '''
    r = re.compile('TCGA-[a-zA-Z\d]{2,4}')
    m = r.search(file)
    if m:
        cancer_type=(str(m.group(0)))
    return cancer_type

def labeler(df,cancer_type):
    labels = []
    index = df.index
    sample_list = list(index)
    for sample in sample_list:
        # extract the last 3 chars in the sample ID column
        r = re.compile('(?<=-)\d\dA')
        m = r.search(str(sample))
        if m:
            if "01" in str(m.group(0)) or "11" in str(m.group(0)):
                labels.append("Tumor")
            else:
                labels.append("Normal")
        else:
            labels.append("remove")
    # Now the labels column to df as the labels column
    df['cancerType'] = cancer_type
    # add the cancertype column
    df['labels'] = labels
    return df

# iterate in 5k blocks, to length of the files
import glob
def build_input(path, start, end):
    block = []
    for file in glob.glob(path):
        dat = read_cut(file, start, end)
        cancer_type = CT(file)
        df = labeler(dat,cancer_type)
        block.append(df)
    df_block = pd.concat(block, ignore_index=True)
    df_block = df_block[df_block.labels != "remove"]

    return df_block

def rfR(X):
    y = X["cancerType"]
    sampleIDs = X["Composite Element REF"]
    X = X.drop("cancerType", axis =1)
    X = X.drop("Composite Element REF", axis =1)
    X = pd.get_dummies(X,columns = ["labels"])
    #X = np.array(X)
    y = np.array(y)
    y = y.reshape(-1, 1)
    encoder = OneHotEncoder(sparse=False)
    encoder.fit(y)
    y_encoded = encoder.transform(y)
    # List of features for later use
    feature_list = list(X.columns)
    X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.3, random_state=19) # 19 love you bubba
    rf = RandomForestRegressor(n_estimators = 100, random_state = 19, n_jobs=20)
    rf.fit(X_train, y_train)
    # Get numerical feature importances
    importances = list(rf.feature_importances_)
    # List of tuples with variable and importance
    feature_importances = [(feature, round(importance, 7)) for feature, importance in zip(feature_list, importances)]
    feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
    df = pd.DataFrame (data = feature_importances[1:50], columns= ["cgIDs","importance_scores"])
    return df
'''
def grade_all(path):
    windows = chunkIt(range(485000), 100)
    for i in windows:
        start = list(i)[0]
        end = list(i)[-1]
        df_block = build_input(path,start,end)
        df = rfR(df_block)
        f = open("window_"+ str(list(windows)[0][0]) + "_" + str(list(windows)[1][0]) + ".csv",'w')
        for index, row in df.iterrows():
            f.write(f'Index: {index}, row: {row.values}')
        f.close()
    pass


grade_all(path)
'''

import time
tic = time.perf_counter()
df_block = build_input(path, 0,5000)
score = rfR(df_block)
toc = time.perf_counter()
print(f"Scored features in {toc - tic:0.4f} seconds")
#windows = chunkIt(range(485000), 100)

