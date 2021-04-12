#!/usr/bin/env python
# coding: utf-8

# import libraries
import io
import math
import os
import re
import time
import warnings

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
import glob

# data loader module
'''
This module will accept normalized beta values input the seelcted features into the trained model
We will then predict the classes of the data. We will report the outcomes as a barchart or pie chart with class confidences
'''

# read in the file with the selected features and make them an object to load:

selected_features = pd.read_csv("Selected_features_all_methylation_data_TCGA.csv")

selected_features.shape


# In[8]:


def feature_IDs(selected_features):
    if "Unnamed: 0" in selected_features.columns:
        selected_features = selected_features.drop("Unnamed: 0", axis=1)
        
    if "labels" in selected_features.columns:
        labels = selected_features.labels
        selected_features = selected_features.drop("labels", axis=1)
    
    if "cancerType" in selected_features.columns:
        cancerType = selected_features.cancerType
        selected_features = selected_features.drop("cancerType", axis=1)

    IDs = selected_features.columns

    return IDs, labels, cancerType

IDs, labels, cancerType = feature_IDs(selected_features)


# In[9]:


IDs


# In[93]:


# read in the data for the normalized beta values and validate this model works on an outside data set
# COADREAD_test
def data_loader(IDs):
    df = pd.read_csv("/home/jovyan/storage/methylation_illumina450k/data/GSE116298_Normalized_betas_GBM.csv")
    
    if "Unnamed: 0" in df.columns:
        df = df[df["Unnamed: 0"].isin(IDs)]
        feature_names = df["Unnamed: 0"]
        df = df.drop("Unnamed: 0", axis=1).T
    # transpose dataframe
    
    df.columns = feature_names
    #df = df.drop("Unnamed: 0", axis=1)
    #df.labels = "Tumor"
    #df.cancerType = "GBM"
    
    
    return df, feature_names

df, feature_names = data_loader(IDs)


# In[94]:


main_list = list(set(selected_features) - set(feature_names))
# wellllll this is soething we have have to re-adjust
for col in main_list:
    if col != "Unnamed: 0":
        df[col] = 0
    
df.labels = "Tumor"
df = pd.get_dummies(df, columns=["labels"])
df.cancerType = "GBM"
df["labels_Normal"] = 0


# In[95]:


df.shape


# In[ ]:


'''# importe required libraries
import openpyxl
import csv
import pandas as pd
  
# open given workbook 
# and store in excel object 
excel = openpyxl.load_workbook("/home/jovyan/storage/methylation_illumina450k/data/GSE166212_Colorectal_Cancer_DNA_Methylation_normalized_betas.xlsx")
  
# select the active sheet
sheet = excel.active
  
# writer object is created
col = csv.writer(open("/home/jovyan/storage/methylation_illumina450k/data/GSE166212_Colorectal_Cancer_DNA_Methylation_normalized_betas.csv",
                      'w', 
                      newline=""))
  
# writing the data in csv file
for r in sheet.rows:
    # row by row write 
    # operation is perform
    col.writerow([cell.value for cell in r])
  
# read the csv file and 
# covert into dataframe object 
#df = pd.DataFrame(pd.read_csv("tt.csv"))
'''


# In[100]:


df = df.drop("cancerType", axis=1)


# In[98]:


# read in the df into the first layer of a NN and begin the grind

import collections
import random

import pandas as pd

# import matplotlib
import tensorflow as tf

print("Num GPUs Available: ", len(tf.config.list_physical_devices("GPU")))
import os

import imblearn
from sklearn.model_selection import cross_val_score
from tensorflow import keras
from tensorflow.keras import Sequential, backend
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier

path = os.getcwd()
print(path)

model = tf.keras.models.load_model('/home/jovyan/storage/Machine_Learning/Illumin450kANNPretrained.h5')


# In[101]:


external_pred = model.predict(df)
y_pred = np.argmax(external_pred, axis=1)
#rounded_labels = np.argmax(y_test, axis=1)
#rounded_labels


# In[111]:


rounded_labels =np.full(
  shape=47,
  fill_value=8,
  dtype=np.int
)


# In[112]:


import numpy as np
import sklearn.metrics as metrics
print("NeuralNet Accuracy: " + str(metrics.accuracy_score(rounded_labels, y_pred)))


# In[ ]:




