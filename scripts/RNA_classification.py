#!/usr/bin/env python
# coding: utf-8

# In[3]:


# lets play with RNA today. and add DNA and RNA to the equation with classes
# pass the classes into the same classifier mathods and get the data back. The we will make mixed classes that accepts incomplete 
# data as input and we will see if we can simulate passing different cases with different information to the total model
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

warnings.filterwarnings("ignore")


# In[55]:


""" Make the reader for feature subsetting

Read in each dataframe, subset into the current row block 5000.
transform the DF to have features across the columns with all patients as row
Add a label to the edge. save file. 
"""

path = "/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/*.csv"
#file = "/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/TCGA-SARC.counts.csv"


def chunkIt(seq, num):
    """
    Gonna split the feature sets into chunks
    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last) : int(last + avg)])
        last += avg

    return out

def read_cut(file, start,end):
    sampleIDs = pd.read_table(file, nrows=0, delimiter =",").columns
    df = (
        pd.read_table(
            file,
            delimiter = ",",
            skiprows=start,
            nrows=end - start,
            names=sampleIDs,
        )
        .fillna(0)
        .set_index("Ensembl")
        .T)
    if "__no_feature" in df.columns:
        df = df.drop("__no_feature", axis=1)
    if "__ambiguous" in df.columns:
        df = df.drop("__ambiguous", axis=1)
    if "__too_Low_aQual" in df.columns:
        df = df.drop("__too_low_aQual", axis=1)
    if "__not_aligned" in df.columns:
        df = df.drop("__not_aligned", axis=1)
    if "__alignment_not_unique" in df.columns:
        df = df.drop("__alignment_not_unique", axis=1)
    
    df.rename(columns = {'Ensembl':'samples'}, inplace = True)
    
    return df

def CT(file):
    """
    Str extract the Cancer type for the
    """
    r = re.compile("TCGA-[a-zA-Z\d]{2,6}")
    m = r.search(file)
    if m:
        cancer_type = str(m.group(0))
    else:
        r = re.compile("TARGET-[a-zA-Z\d]{2,6}")
        m = r.search(file)
        cancer_type = str(m.group(0))
    return cancer_type

def labeler(df, cancer_type):
    labels = []
    index = df.index
    sample_list = list(index)
    for sample in sample_list:
        # extract the last 3 chars in the sample ID column
        r = re.compile("(?<=-)\d\d")
        m = r.search(str(sample))
        if m:
            if "01" in str(m.group(0)) or "11" in str(m.group(0)):
                labels.append("Tumor")
            else:
                labels.append("Normal")
        else:
            labels.append("remove")
    # Now the labels column to df as the labels column
    df["cancerType"] = cancer_type
    # add the cancertype column
    df["labels"] = labels
    return df

def build_input(path, start, end):
    block = []
    for file in glob.glob(path):
        dat = read_cut(file, start, end)
        cancer_type = CT(file)
        df = labeler(dat, cancer_type)
        block.append(df)
    df_block = pd.concat(block, ignore_index=True)
    df_block = df_block[df_block.labels != "remove"]

    return df_block

def rfR(X):
    y = X["cancerType"]
    #sampleIDs = X["samples"]
    X = X.drop("cancerType", axis=1)
    if "Composite Element REF" in X.columns:
        X = X.drop("Composite Element REF", axis=1)
    if "samples" in X.columns:
        X = X.drop("samples", axis=1)
    X = pd.get_dummies(X, columns=["labels"])
    y = np.array(y)
    y = y.reshape(-1, 1)
    encoder = OneHotEncoder(sparse=False)
    encoder.fit(y)
    y_encoded = encoder.transform(y)
    # List of features for later use
    feature_list = list(X.columns)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_encoded, test_size=0.3, random_state=19
    )  # 19 love you bubba
    rf = RandomForestRegressor(n_estimators=100, random_state=19, n_jobs=20)
    rf.fit(X_train, y_train)
    # Get numerical feature importances
    importances = list(rf.feature_importances_)
    # List of tuples with variable and importance
    feature_importances = [
        (feature, round(importance, 7))
        for feature, importance in zip(feature_list, importances)
    ]
    feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
    df = pd.DataFrame(
        data=feature_importances[1:50], columns=["Ensembl", "importance_scores"]
    )
    return df


def grade_all(path):
    windows = chunkIt(range(60486), 100)
    for i in windows:
        start = list(i)[0]
        end = list(i)[-1]
        print(start)
        print(end)
        df_block = build_input(path, start, end)
        score = rfR(df_block)
        score.to_csv("Ensembl_ID_window_" + str(start) + "_" + str(end) + ".csv")
    pass


# In[42]:


start = 59881
end = 60486
dat = read_cut(file, start =59881, end = 60486)
cancer_type = CT(file)
df = labeler(dat, cancer_type)
df_block = build_input(path, start =59881, end = 60486)
score = rfR(df_block)
score.to_csv("Ensembl_ID_window_" + str(start) + "_" + str(end) + ".csv")


# In[37]:


windows = chunkIt(range(60486), 100)


# In[38]:


df =grade_all(path)


# In[32]:


pwd


# In[47]:


path = "/home/jovyan/storage/Machine_Learning/RNA_ID_windows/*.csv"


# In[50]:


def join(path):
    features = []
    for file in glob.glob(path):
        dat = pd.read_csv(file)
        features.append(dat)
        # concatnenate them and sort the feature ranks
    df = pd.concat(features, ignore_index=True)
    df.drop_duplicates(subset=["Ensembl"])
    df.sort_values(by="importance_scores", ascending=False)
    # get cgIDs
    features = df.Ensembl
    features = features.drop_duplicates()
    return features


# In[57]:


features


# In[ ]:


import time

def subset(path, features):
    selected = []
    tic_all = time.perf_counter()
    for file in glob.glob(path):
        print(file)
        df = pd.read_table(file, delimiter=",")
        df = df[df["Ensembl"].isin(features)]
        df = df.T
        names = df.iloc[0]
        cancer_type = CT(file)
        df = labeler(df, cancer_type)
        selected.append(df)
    df_block = pd.concat(selected, ignore_index=True)
    df_block = df_block.fillna(0)
    df_block = df_block[df_block.labels != "remove"]
    print(df_block.shape)
    column_indices = range(4702)
    new_names = names
    old_names = df_block.columns[column_indices]
    df_block.rename(columns=dict(zip(old_names, new_names)), inplace=True)
    toc_all = time.perf_counter()
    print(f"Read in features in {toc_all - tic_all:0.4f} seconds")
    return df_block, names

features = join(path)
df, names = subset(path, features)
print(df.columns[1:10])
print(names[1:10])


# In[ ]:


df


# In[ ]:


df.to_csv("Selected_features_all_RNA_counts_TCGA_2.csv")


# In[ ]:




