# import libraries for the gradient boosting
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xgboost
import math
from __future__ import division
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import explained_variance_score
from xgboost import XGBRegressor
from matplotlib import pyplot
import re

# Make the reader for feature subsetting
'''
Read in each dataframe, subset into the current row block 5000.
transform the DF to have features across the columns with all patients as row
Add a label to the edge. save file. 
'''

path = os.listdir("/home/jovyan/storage/CUP/tsv/Compressed/")


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



def reader(file):
    '''
    Read in Methylation data
    '''
    dat = pd.read_csv(file,delimiter = "\t", compression="gzip")
    dat = dat.fillna(0)
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


def cut(dat, start, end):
    df = dat[start:end:]
    df = df.T
    cgIDs = df.iloc[0] #grab the first row for the header
    df = df[1:] #take the data less the header row
    df.columns = cgIDs #set the header row as the df header
    return df

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
#############
#dummy runs
#############
dat = reader("~/storage/CUP/tsv/Compressed/TCGA-KIRC.methylation450.tsv.gz")
cancer_type = CT("~/storage/CUP/tsv/Compressed/TCGA-KIRC.methylation450.tsv.gz")
df = cut(dat, 0, 5000)
df = labeler(df,cancer_type)

###########

# build input to score 5000 features
import glob
def build_input(path_list, start,end):
    block = []
    #start = 0
    #end = 5000
    for file in path_list:
        print(str(file))
        dat = reader(file)
        cancer_type = CT(file)
        df = cut(dat, start, end)
        df = labeler(df,cancer_type)
        block.append(df)
    df_block = pd.concat(block, ignore_index=True)
    df_block = df_block[df_block.labels != "remove"]
    return df_block

df_block = build_input(path_list = glob.glob("/home/jovyan/storage/CUP/tsv/Compressed/*.gz"), start = 0, end = 5000)
write.csv(df_block, "/home/jovyan/storage/CUP/challenge_data.csv")

############
# iterate over the window widths and make all of the splits
def grade_all(path_list)
  windows = chunkIt(range(485000), 100)
  for window in windows:
    start = list(window)[0]
    end = list(window)[-1]
    df = build_input(path_list = glob.glob("/home/jovyan/storage/CUP/tsv/Compressed/*.gz"), start, end)
    
    # make model for these 5000
    
    # reduce the data set to top 1% 
    
    # Add to final list. 
        
        
        
        
        
        
